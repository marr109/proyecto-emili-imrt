#!/usr/bin/env bash
set -uo pipefail

# Uso: scripts/run_sa_budget_pilot.sh [segundos] [paralelismo] [T_inicial] [T_final]
#   por defecto: 1080 s (18 min), 4 simultaneas, T de 1300 a 80
#
# MUESTRA PILOTO de recocido simulado sobre BAO con evaluacion exacta del FMO.
#
# Que responde y que NO responde
# ------------------------------
# Responde: si SA con enfriamiento es funcional sobre este problema y donde
# queda su objetivo final frente a ILS bajo el mismo presupuesto de 18 minutos.
#
# NO responde cual es el mejor esquema de enfriamiento. Barrer temperaturas a
# mano seria hacer a mano lo que irace hace sistematicamente, y el informe ya
# declara a irace como el mecanismo de calibracion.
#
# TRAMPA DEL PARSER: "ratio" no es una razon, es un decremento absoluto
# ---------------------------------------------------------------------
# emili::Metropolis::accept (emilibase.cpp) actualiza la temperatura asi:
#
#     temperature = (alpha * temperature) - beta;     // beta == el token "ratio"
#
# El constructor de tres argumentos que usa sa_metropolis fija alpha = 1 y
# interval = 1, de modo que el esquema resultante es LINEAL y se aplica en cada
# iteracion: T_{k+1} = T_k - beta. Pasarle 0.95 esperando enfriamiento
# geometrico restaria 0.95 grados por iteracion y la temperatura no bajaria
# nunca dentro del presupuesto. Para geometrico real hace falta el token saacc
# con ratio 0 y alpha < 1.
#
# De donde salen T_inicial, T_final y el decremento
# -------------------------------------------------
# Metropolis acepta un empeoramiento con probabilidad exp(-delta/T), asi que T
# solo tiene sentido en la escala de los delta del problema. Midiendo los 7387
# empeoramientos observados en las 120 trayectorias de la matriz de presupuesto
# (experiments/ils/*/**-budget1080s/):
#
#     p25 = 360    mediana = 901    p75 = 1919    media = 1427
#
# T_inicial = 1300: exp(-901/1300) = 0,50. Al arrancar acepta la mitad de los
# empeoramientos de magnitud tipica, que es el regimen de exploracion.
#
# T_final = 80: exp(-360/80) = 0,011. Al cerrar rechaza practicamente todo,
# incluidos los empeoramientos chicos, y la busqueda queda en descenso puro.
#
# El decremento sale del presupuesto REAL de iteraciones, no de una convencion,
# porque el framework solo sabe enfriar por iteracion y no por tiempo.
#
# La matriz de ILS rendia 14 evaluaciones por minuto, pero ese numero NO aplica
# aca: buena parte de esas evaluaciones pegaban en la cache del vecindario y no
# llamaban a Gurobi. SA muestrea con Neighborhood::random, falla la cache casi
# siempre y paga el solve completo en cada iteracion. Medido directamente sobre
# una corrida SA de 120 s (trajectory.csv, columna cached en false en todas las
# filas): 12 evaluaciones, es decir 6 por minuto, unas 108 en los 18 minutos.
#
# Esa medicion se tomo con otros cuatro procesos compitiendo por la maquina y
# con threads=2, asi que 6 por minuto es un PISO y no el rendimiento esperado
# del piloto, que correra solo. Aun asi se calibra contra el piso, por la
# asimetria del error:
#
#   - Si se calibra al techo y el rendimiento real es el piso, la temperatura
#     se queda en ~760 al agotarse el presupuesto: nunca enfria, y el algoritmo
#     degenera en la caminata aleatoria sesgada que este cambio vino a eliminar.
#   - Si se calibra al piso y el rendimiento real es el techo, la temperatura
#     llega a T_final antes de tiempo y el resto de la corrida es descenso puro,
#     que es justamente el cierre que se le pide a un esquema de enfriamiento.
#
# Enfriar de 1300 a 80 en 108 pasos pide beta = (1300 - 80) / 108 = 11,3.

BUDGET="${1:-1080}"
PAR="${2:-4}"
T_START="${3:-1300}"
T_END="${4:-80}"
STEP="${STEP:-5}"
# NEIGH selecciona el operador que SA usa para muestrear:
#   shift  -> nangshift <STEP>. random() sortea saltos en (STEP, 2*STEP): con
#             STEP=5 son 6..9 posiciones. Es el operador de PERTURBACION de ILS,
#             no el vecindario que enumeran First/Best. Rompe la invariante de
#             residuo modulo STEP, que es lo que lo hace ergodico, pero su salto
#             minimo es mayor que un paso de busqueda local: al enfriar no puede
#             refinar porque todos sus candidatos caen lejos.
#   multi  -> nangshiftmulti con escala 1 2 3 5 8 13 21. Sortea el paso de la
#             lista en vez de inflarlo, asi que cubre grueso y fino con el mismo
#             operador. Es ergodico porque la lista incluye el 1. La seleccion
#             entre escalas la hace la propia aceptacion: caliente admite saltos
#             de 21, frio solo sobrevive el paso de 1.
NEIGH="${NEIGH:-shift}"
CATALOG="${CATALOG:-restricted}"
SEEDS="${SEEDS:-15}"
# Iteraciones esperadas dentro del presupuesto, al piso medido de 6 por minuto.
ITERS="${ITERS:-$(( BUDGET * 6 / 60 ))}"
BETA="${BETA:-$(awk -v a="$T_START" -v b="$T_END" -v n="$ITERS" \
  'BEGIN{ v=(a-b)/n; printf (v<1 ? "%.3f" : "%.0f"), v }')}"

if [ -z "${EMILI_AMPL_BIN_DIR:-}" ] || [ -z "${EMILI_GUROBI_BIN:-}" ]; then
  SP="$(echo ampl_gurobi/.venv/lib/python*/site-packages)"
  [ -d "$SP/ampl_module_base" ] || { echo "No se encontro ampl_module_base bajo $SP" >&2; exit 1; }
  export EMILI_AMPL_BIN_DIR="$PWD/$SP/ampl_module_base/bin"
  export EMILI_GUROBI_BIN="$PWD/$SP/ampl_module_gurobi/bin/gurobi"
fi
export EMILI_GUROBI_OPTIONS="${EMILI_GUROBI_OPTIONS:-threads=4}"

case "$NEIGH" in
  shift) NEIGH_TOKEN="nangshift $STEP"; NEIGH_SUF="" ;;
  multi) NEIGH_TOKEN="nangshiftmulti 7 1 2 3 5 8 13 21"; NEIGH_SUF="-multi" ;;
  *) echo "NEIGH debe ser shift o multi (recibido: $NEIGH)" >&2; exit 1 ;;
esac

PY="ampl_gurobi/.venv/bin/python3"
CAT_TOKEN=""; [ "$CATALOG" = "unrestricted" ] && CAT_TOKEN="unrestricted"
OUT="experiments/sa/nangshift${STEP}${NEIGH_SUF}/${CATALOG}-budget${BUDGET}s-T${T_START}to${T_END}/data"
MASTER="experiments/sa/pilot${NEIGH_SUF}_T${T_START}to${T_END}_budget${BUDGET}s.log"
mkdir -p "$OUT" "$(dirname "$MASTER")"

log() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a "$MASTER"; }

run_one() {
  local seed="$1"
  local dir="${OUT}/seed${seed}"
  local traj="${dir}/trajectory.csv" rlog="${dir}/run.log"
  mkdir -p "$dir"

  # El corte por presupuesto sale por finalise()/exit(0) y saltea el bloque de
  # main() que imprime "Found solution"; el reporte clinico (atexit) es lo unico
  # presente en los dos caminos.
  if [ -f "$rlog" ] && grep -q "Conformity Index" "$rlog" 2>/dev/null; then
    return 0
  fi

  # tmaxiter inalcanzable: el unico criterio de corte es -it, igual que en la
  # matriz de ILS, para que ambos compitan por el mismo recurso.
  # -it va ANTES de rnds: getTime() usa checkToken() sobre la posicion actual
  # del cursor mientras rnds usa seek(), asi que al final se ignora en silencio.
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    sa irandomk $CAT_TOKEN \
    tmaxiter 1000000 \
    $NEIGH_TOKEN \
    sa_metropolis "$T_START" "$T_END" "$BETA" \
    -it "$BUDGET" \
    rnds "$seed" \
    > "$rlog" 2>&1

  if grep -q "Conformity Index" "$rlog" 2>/dev/null; then
    "$PY" scripts/plot_trajectory.py --csv "$traj" --skip-trajectory \
      --improvements-csv "${dir}/trajectory_improvements.csv" \
      --convergence-output "${dir}/convergence.png" \
      --iterations-output "${dir}/iterations.png" >> "$rlog" 2>&1 || true
  else
    echo "WARNING seed${seed}: sin marca de corrida completa" >> "$MASTER"
  fi
}
export -f run_one
export BUDGET T_START T_END BETA STEP NEIGH NEIGH_TOKEN NEIGH_SUF CATALOG CAT_TOKEN OUT PY MASTER \
       EMILI_AMPL_BIN_DIR EMILI_GUROBI_BIN EMILI_GUROBI_OPTIONS

log "PILOTO SA START  vecindario=${NEIGH} [${NEIGH_TOKEN}]  T=${T_START}->${T_END} beta=${BETA} (${ITERS} iter estimadas)  paso=${STEP}  catalogo=${CATALOG}  presupuesto=${BUDGET}s  semillas=${SEEDS}  paralelismo=${PAR}"

DONE=0
for seed in $(seq -w 1 "$SEEDS"); do
  run_one "$seed" &
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
  DONE=$((DONE + 1))
  [ $((DONE % PAR)) -eq 0 ] && log "lanzadas ${DONE}/${SEEDS}"
done
wait

COMP=$(find "$OUT" -name run.log -exec grep -l "Conformity Index" {} \; 2>/dev/null | wc -l)
log "PILOTO SA COMPLETE  corridas completas: ${COMP}/${SEEDS}"
