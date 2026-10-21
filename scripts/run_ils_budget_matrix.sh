#!/usr/bin/env bash
set -uo pipefail

# Uso: scripts/run_ils_budget_matrix.sh [segundos] [paralelismo]
#   por defecto: 1080 s (18 min) y 4 corridas simultaneas
#
# Variables de entorno:
#   ORDERS      ordenes de escaneo a correr: "fijo random circular" (default: fijo)
#   STRATEGIES  estrategias: "first best" (default: ambas)
#   SEEDS       cantidad de semillas (default: 15)
#
# Corre la matriz completa de ILS bajo presupuesto de TIEMPO:
#   {restricted, unrestricted} x {nangshift 5, 10} x {First, Best} x 15 semillas
#   = 8 condiciones, 120 corridas.
#
# Salida:
#   experiments/ils/nangshift<step>/<catalog>-budget<secs>s/data/<strategy>/seedNN/
#
# Decisiones de diseño, medidas antes de fijarlas:
#
#   1. tmaxiter va a un valor inalcanzable y el corte lo da SOLO -it. Si ambos
#      criterios pudieran cortar, First podria completar sus iteraciones sin
#      agotar el reloj mientras a Best lo corta el timer: ni igual tiempo ni
#      iguales iteraciones. El presupuesto de tiempo debe ser la unica variable.
#
#   2. -it va ANTES de rnds. getTime() usa tm.checkToken(), que mira el token en
#      la posicion actual del cursor; rnds usa tm.seek(), que busca en cualquier
#      lado. Un -it despues de rnds se ignora EN SILENCIO.
#
#   3. Paralelismo 4. Medido en calafate ociosa con presupuesto de 120 s:
#        serial  -> 15 solves
#        4 par.  -> 14 14 14 14      (uniforme)
#        8 par.  -> 11..13           (dispersion 18%)
#        16 par. -> 10..12
#      La contencion no es de CPU: bajar los threads de Gurobi de 4 a 1 no la
#      mueve. Es ancho de banda de memoria por el marshalling de ~850k tuplas
#      por solve. Con 4 el presupuesto efectivo es identico entre corridas, que
#      es lo que hace comparables las condiciones.

BUDGET="${1:-1080}"
PAR="${2:-4}"
SEEDS="${SEEDS:-15}"
ORDERS="${ORDERS:-fijo}"
STRATEGIES="${STRATEGIES:-first best}"

# El orden de escaneo solo altera el resultado en First: Best recorre el
# vecindario completo antes de moverse, asi que el orden unicamente puede
# resolver empates. El token va pegado a nangshift <step>.
orden_token() { case "$1" in random) echo "randorder";; circular) echo "circular";; *) echo "";; esac; }
orden_sufijo() { case "$1" in random) echo "-rand";; circular) echo "-circ";; *) echo "";; esac; }

if [ -z "${EMILI_AMPL_BIN_DIR:-}" ] || [ -z "${EMILI_GUROBI_BIN:-}" ]; then
  SP="$(echo ampl_gurobi/.venv/lib/python*/site-packages)"
  [ -d "$SP/ampl_module_base" ] || { echo "No se encontro ampl_module_base bajo $SP" >&2; exit 1; }
  export EMILI_AMPL_BIN_DIR="$PWD/$SP/ampl_module_base/bin"
  export EMILI_GUROBI_BIN="$PWD/$SP/ampl_module_gurobi/bin/gurobi"
fi
export EMILI_GUROBI_OPTIONS="${EMILI_GUROBI_OPTIONS:-threads=4}"

PY="ampl_gurobi/.venv/bin/python3"
MASTER="experiments/ils/budget${BUDGET}s_matrix.log"
mkdir -p "$(dirname "$MASTER")"

log() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$MASTER"; }

run_one() {
  local step="$1" cat="$2" strat="$3" seed="$4" orden="${5:-fijo}"
  local tok=""; [ "$cat" = "unrestricted" ] && tok="unrestricted"
  local otok="$(orden_token "$orden")"
  local osuf="$(orden_sufijo "$orden")"
  local dir="experiments/ils/nangshift${step}/${cat}-budget${BUDGET}s${osuf}/data/${strat}/seed${seed}"
  local traj="${dir}/trajectory.csv" rlog="${dir}/run.log"
  mkdir -p "$dir"

  # El corte por presupuesto sale via finalise()/exit(0), que saltea el bloque
  # de main() donde se imprime "Found solution". El reporte clinico (atexit) es
  # lo unico presente en los dos caminos, asi que esa es la marca de completa.
  if [ -f "$rlog" ] && grep -q "Conformity Index" "$rlog" 2>/dev/null; then
    return 0
  fi

  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    ils "$strat" irandomk $tok locmin nangshift "$step" $otok \
    tmaxiter 1000000 \
    prangshift "$step" 3 \
    baoimprove rejectrepeated \
    -it "$BUDGET" \
    rnds "$seed" \
    > "$rlog" 2>&1

  if grep -q "Conformity Index" "$rlog" 2>/dev/null; then
    "$PY" scripts/plot_trajectory.py --csv "$traj" --skip-trajectory \
      --improvements-csv "${dir}/trajectory_improvements.csv" \
      --convergence-output "${dir}/convergence.png" \
      --iterations-output "${dir}/iterations.png" >> "$rlog" 2>&1 || true
  else
    echo "WARNING ${step}/${cat}${osuf}/${strat}/seed${seed}: sin marca de corrida completa" >> "$MASTER"
  fi
}
export -f run_one orden_token orden_sufijo
export BUDGET PY MASTER ORDERS STRATEGIES EMILI_AMPL_BIN_DIR EMILI_GUROBI_BIN EMILI_GUROBI_OPTIONS

# --- cola de trabajos -------------------------------------------------------
JOBS="$(mktemp)"
for step in 5 10; do
  for cat in restricted unrestricted; do
    for orden in $ORDERS; do
      for strat in $STRATEGIES; do
        for seed in $(seq -w 1 "$SEEDS"); do
          echo "$step $cat $strat $seed $orden" >> "$JOBS"
        done
      done
    done
  done
done
TOTAL=$(wc -l < "$JOBS")

log "MATRIZ START  presupuesto=${BUDGET}s  paralelismo=${PAR}  corridas=${TOTAL}"
log "ordenes=[${ORDERS}]  estrategias=[${STRATEGIES}]  semillas=${SEEDS}"
log "gurobi_options=${EMILI_GUROBI_OPTIONS}"
log "estimado: $(( (TOTAL / PAR) * (BUDGET + 30) / 3600 )) h aprox"

DONE=0
while read -r step cat strat seed orden; do
  run_one "$step" "$cat" "$strat" "$seed" "$orden" &
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
  DONE=$((DONE + 1))
  [ $((DONE % PAR)) -eq 0 ] && log "lanzadas ${DONE}/${TOTAL}"
done < "$JOBS"
wait
rm -f "$JOBS"

log "MATRIZ COMPLETE"
COMP=$(find experiments/ils -path "*budget${BUDGET}s*" -name run.log -exec grep -l "Conformity Index" {} \; 2>/dev/null | wc -l)
log "corridas completas: ${COMP}/${TOTAL}"
