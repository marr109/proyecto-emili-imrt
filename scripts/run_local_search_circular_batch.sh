#!/usr/bin/env bash
set -euo pipefail

# Uso: scripts/run_local_search_circular_batch.sh <step> <strategy: first|best> <catalog: restricted|unrestricted>
# Corre 15 semillas de búsqueda local pura (sin ILS) con el vecindario nangshift <step> en
# modo circular (arranca cada ronda en el slot recién modificado, rota circularmente),
# contra el catálogo indicado. Salida en:
#   experiments/local_search/nangshift<step>/<catalog>-shift-circular/data/<strategy>/seedNN/
#
# NO toca experiments/local_search/nangshift<step>/restricted/, .../unrestricted/
# ni .../<catalog>-shift-rand/ (baselines fijo y random, no se sobreescriben).

STEP="$1"
STRATEGY="$2"
CATALOG="$3"
PY="ampl_gurobi/.venv/bin/python3"

if [ "$CATALOG" = "unrestricted" ]; then
  CATALOG_TOKEN="unrestricted"
elif [ "$CATALOG" = "restricted" ]; then
  CATALOG_TOKEN=""
else
  echo "catalog debe ser 'restricted' o 'unrestricted', recibido: $CATALOG" >&2
  exit 1
fi

OUTBASE="experiments/local_search/nangshift${STEP}/${CATALOG}-shift-circular/data/${STRATEGY}"
LOG="experiments/local_search/nangshift${STEP}/${CATALOG}-shift-circular/${STRATEGY}_batch.log"

mkdir -p "$OUTBASE"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"

  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya completa: Found solution en run.log)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG})" | tee -a "$LOG"

    ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
      "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" circular \
      rnds "$seed" \
      > "$run_log" 2>&1

    if grep -q "Found solution" "$run_log" 2>/dev/null; then
      echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
    else
      echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no se encontró 'Found solution' en run.log, revisar" | tee -a "$LOG"
      continue
    fi
  fi

  if [ -f "${seed_dir}/convergence.png" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP plot seed${seed} (ya existe convergence.png)" | tee -a "$LOG"
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

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG})" | tee -a "$LOG"
