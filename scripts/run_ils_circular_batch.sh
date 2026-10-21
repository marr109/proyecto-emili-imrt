#!/usr/bin/env bash
set -euo pipefail

# Uso: scripts/run_ils_circular_batch.sh <step> <first|best> <restricted|unrestricted> [tmaxiter]
# Igual que run_ils_batch.sh, pero con el vecindario en círculo
# (nangshift <step> circular) en vez de orden fijo. Salida en:
#   experiments/ils/nangshift<step>/<catalog>-circular/data/<strategy>/seedNN/
#
# tmaxiter es opcional, default 10 (igual que la línea base fija). Si se pasa
# un valor distinto, queda reflejado en el nombre de carpeta para no mezclar
# presupuestos distintos bajo el mismo path.
#
# NO toca las carpetas de orden fijo (grid5/, restricted/, unrestricted/).

STEP="$1"
STRATEGY="$2"
CATALOG="$3"
TMAXITER="${4:-10}"
PY="ampl_gurobi/.venv/bin/python3"

if [ "$CATALOG" = "unrestricted" ]; then
  CATALOG_TOKEN="unrestricted"
elif [ "$CATALOG" = "restricted" ]; then
  CATALOG_TOKEN=""
else
  echo "catalog debe ser 'restricted' o 'unrestricted', recibido: $CATALOG" >&2
  exit 1
fi

if [ "$TMAXITER" = "10" ]; then
  COND_SUFFIX="${CATALOG}-circular"
else
  COND_SUFFIX="${CATALOG}-circular-tmax${TMAXITER}"
fi
OUTBASE="experiments/ils/nangshift${STEP}/${COND_SUFFIX}/data/${STRATEGY}"
LOG="experiments/ils/nangshift${STEP}/${COND_SUFFIX}/${STRATEGY}_batch.log"

mkdir -p "$OUTBASE"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"

  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya completa)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} circular)" | tee -a "$LOG"

    ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
      ils "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" circular \
      tmaxiter "$TMAXITER" \
      prangshift "$STEP" 3 \
      baoimprove rejectrepeated \
      rnds "$seed" \
      > "$run_log" 2>&1

    if grep -q "Found solution" "$run_log" 2>/dev/null; then
      echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
    else
      echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: sin 'Found solution', revisar" | tee -a "$LOG"
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

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} circular)" | tee -a "$LOG"
