#!/usr/bin/env bash
set -euo pipefail

NANGSHIFT=5
OUTBASE="experiments/ils/perturbation_magnitude_pilot/data"
LOG="experiments/ils/perturbation_magnitude_pilot/pilot_batch.log"
PY="ampl_gurobi/.venv/bin/python3"

mkdir -p "$OUTBASE"

for mag in 15 30 45; do
  for seed in 01 02 03; do
    seed_dir="${OUTBASE}/mag${mag}/seed${seed}"
    traj="${seed_dir}/trajectory.csv"

    if [ -f "$traj" ] && [ "$(grep -c "^" "$traj" 2>/dev/null)" -gt 1 ] && [ -f "${seed_dir}/convergence.png" ]; then
      echo "[$(date '+%H:%M:%S')] SKIP mag${mag}/seed${seed} (ya existe)" | tee -a "$LOG"
      continue
    fi

    mkdir -p "$seed_dir"

    if [ -f "$traj" ] && [ "$(grep -c "^" "$traj" 2>/dev/null)" -gt 1 ]; then
      echo "[$(date '+%H:%M:%S')] REPLOT mag${mag}/seed${seed} (trajectory ya existe, faltan plots)" | tee -a "$LOG"
    else
      echo "[$(date '+%H:%M:%S')] RUN mag${mag}/seed${seed}" | tee -a "$LOG"

      ./build/emili instances/CERR_Prostate baoimrt 4 csv "${traj}" \
        ils first irandomk locmin nangshift "$NANGSHIFT" \
        tmaxiter 10 \
        prangshift "$mag" 3 \
        baoimprove rejectrepeated \
        rnds "$seed" \
        >> "$LOG" 2>&1
    fi

    "$PY" scripts/plot_trajectory.py \
      --csv "${traj}" \
      --skip-trajectory \
      --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
      --convergence-output "${seed_dir}/convergence.png" \
      --iterations-output "${seed_dir}/iterations.png" \
      >> "$LOG" 2>&1

    "$PY" scripts/plot_ils_trajectory.py \
      --csv "${traj}" \
      --step "$NANGSHIFT" \
      --output "${seed_dir}/ils_perturbations.png" \
      --perturbations-csv "${seed_dir}/perturbations.csv" \
      --moves-csv "${seed_dir}/perturbation_moves.csv" \
      --title "mag${mag} seed${seed}" \
      >> "$LOG" 2>&1

    echo "[$(date '+%H:%M:%S')] DONE mag${mag}/seed${seed}" | tee -a "$LOG"
  done
done

echo "[$(date '+%H:%M:%S')] PILOT COMPLETE" | tee -a "$LOG"
