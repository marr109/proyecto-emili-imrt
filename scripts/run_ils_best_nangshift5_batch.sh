#!/usr/bin/env bash
set -euo pipefail

STEP=5
COND="nangshift${STEP}/grid5"
OUTBASE="experiments/ils/${COND}/data/best"
LOG="experiments/ils/${COND}/data/best_batch.log"
PY="ampl_gurobi/.venv/bin/python3"

mkdir -p "$OUTBASE"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"

  if [ -f "$traj" ] && [ "$(grep -c "^" "$traj" 2>/dev/null)" -gt 1 ] && [ -f "${seed_dir}/convergence.png" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya existe: $traj + plots)" | tee -a "$LOG"
    continue
  fi

  mkdir -p "$seed_dir"

  if [ -f "$traj" ] && [ "$(grep -c "^" "$traj" 2>/dev/null)" -gt 1 ]; then
    echo "[$(date '+%H:%M:%S')] REPLOT seed${seed} (trajectory ya existe, faltan plots)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed}" | tee -a "$LOG"

    ./build/emili instances/CERR_Prostate baoimrt 4 csv "${seed_dir}/trajectory.csv" \
      ils best irandomk locmin nangshift "$STEP" \
      tmaxiter 10 \
      prangshift "$STEP" 3 \
      baoimprove rejectrepeated \
      rnds "$seed" \
      >> "$LOG" 2>&1
  fi

  "$PY" scripts/plot_trajectory.py \
    --csv "${seed_dir}/trajectory.csv" \
    --skip-trajectory \
    --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
    --convergence-output "${seed_dir}/convergence.png" \
    --iterations-output "${seed_dir}/iterations.png" \
    >> "$LOG" 2>&1

  echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
done

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE" | tee -a "$LOG"
