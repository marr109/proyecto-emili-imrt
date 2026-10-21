#!/usr/bin/env bash
set -euo pipefail

# Uso: scripts/run_ils_budget_batch.sh <step> <strategy: first|best> <catalog: restricted|unrestricted> <segundos>
#
# Corre 15 semillas de ILS con PRESUPUESTO DE TIEMPO en vez de un número fijo de
# perturbaciones, para comparar estrategias bajo el mismo esfuerzo de reloj.
# Salida en:
#   experiments/ils/nangshift<step>/<catalog>-budget<segundos>s/data/<strategy>/seedNN/
#
# Diferencias con run_ils_batch.sh, y por qué:
#
#   1. tmaxiter queda en un valor efectivamente inalcanzable y el corte lo da -it.
#      El ILS necesita un objeto de terminación sí o sí (buildTermination devuelve
#      nulo si no matchea ningún token), así que no se puede simplemente omitir.
#
#   2. -it va ANTES de rnds. getTime() usa tm.checkToken(), que mira el token en la
#      posición actual del cursor, mientras que rnds usa tm.seek(), que busca en
#      cualquier lado. Si -it queda después de rnds el parser lo ignora EN SILENCIO
#      y la corrida se va sin presupuesto.
#
#   3. No se llama a plot_ils_trajectory.py. Ese script reconstruye las
#      perturbaciones simulando la semántica del ILS y su propia documentación
#      advierte "Solo válido para `first` (Best no rebasa a mitad de ronda)".
#      Usarlo sobre trayectorias de best daría una clasificación incorrecta.

STEP="$1"
STRATEGY="$2"
CATALOG="$3"
BUDGET="$4"

if [ "$CATALOG" = "unrestricted" ]; then
  CATALOG_TOKEN="unrestricted"
elif [ "$CATALOG" = "restricted" ]; then
  CATALOG_TOKEN=""
else
  echo "catalog debe ser 'restricted' o 'unrestricted', recibido: $CATALOG" >&2
  exit 1
fi

if [ "$STRATEGY" != "first" ] && [ "$STRATEGY" != "best" ]; then
  echo "strategy debe ser 'first' o 'best', recibido: $STRATEGY" >&2
  exit 1
fi

# El fallback compilado en imrt_fmo.cpp apunta a una ruta con python3.9, que no
# existe en todos los entornos. Se resuelve el venv real una sola vez.
if [ -z "${EMILI_AMPL_BIN_DIR:-}" ] || [ -z "${EMILI_GUROBI_BIN:-}" ]; then
  SP="$(echo ampl_gurobi/.venv/lib/python*/site-packages)"
  if [ ! -d "$SP/ampl_module_base" ]; then
    echo "No se encontró ampl_module_base bajo $SP." >&2
    echo "Exportá EMILI_AMPL_BIN_DIR y EMILI_GUROBI_BIN a mano." >&2
    exit 1
  fi
  export EMILI_AMPL_BIN_DIR="$PWD/$SP/ampl_module_base/bin"
  export EMILI_GUROBI_BIN="$PWD/$SP/ampl_module_gurobi/bin/gurobi"
fi

PY="ampl_gurobi/.venv/bin/python3"
OUTBASE="experiments/ils/nangshift${STEP}/${CATALOG}-budget${BUDGET}s/data/${STRATEGY}"
LOG="experiments/ils/nangshift${STEP}/${CATALOG}-budget${BUDGET}s/${STRATEGY}_batch.log"

mkdir -p "$OUTBASE" "$(dirname "$LOG")"

echo "[$(date '+%H:%M:%S')] BATCH START step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} budget=${BUDGET}s" | tee -a "$LOG"
echo "[$(date '+%H:%M:%S')] AMPL bin: ${EMILI_AMPL_BIN_DIR}" | tee -a "$LOG"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"

  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Conformity Index" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya completa)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} budget=${BUDGET}s)" | tee -a "$LOG"

    # tmaxiter 1000000: tope inalcanzable, el corte real lo da -it.
    # -it antes de rnds: obligatorio, ver nota 2 en la cabecera.
    ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
      ils "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" \
      tmaxiter 1000000 \
      prangshift "$STEP" 3 \
      baoimprove rejectrepeated \
      -it "$BUDGET" \
      rnds "$seed" \
      > "$run_log" 2>&1

    # El corte por presupuesto termina via finalise()/exit(0), que NO imprime
    # "Found solution": ese bloque de main() queda saltado. Lo unico que aparece
    # en los dos caminos es el reporte clinico, emitido por atexit, asi que esa
    # es la marca de corrida completa.
    if grep -q "Conformity Index" "$run_log" 2>/dev/null; then
      echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
    else
      echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: sin marca de corrida completa, revisar" | tee -a "$LOG"
      continue
    fi
  fi

  if [ -f "${seed_dir}/convergence.png" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP plot seed${seed}" | tee -a "$LOG"
    continue
  fi

  "$PY" scripts/plot_trajectory.py \
    --csv "$traj" \
    --skip-trajectory \
    --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
    --convergence-output "${seed_dir}/convergence.png" \
    --iterations-output "${seed_dir}/iterations.png" \
    >> "$run_log" 2>&1

  echo "[$(date '+%H:%M:%S')] PLOT seed${seed}" | tee -a "$LOG"
done

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} budget=${BUDGET}s)" | tee -a "$LOG"
