#!/usr/bin/env bash
set -euo pipefail

# Uso: scripts/run_sa_batch.sh
# Corre 15 semillas de recocido simulado (SimulatedAnnealing, aceptación
# Metrópolis) sobre CERR_Prostate, catálogo restricted, nangshift 5,
# tmaxiter 200, sa_metropolis(start=50000, end=1, ratio=250).
# Salida en: experiments/sa/nangshift5/restricted/data/seedNN/

PY="ampl_gurobi/.venv/bin/python3"
OUTBASE="experiments/sa/nangshift5/restricted/data"
LOG="experiments/sa/nangshift5/restricted/batch.log"

mkdir -p "$OUTBASE"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"

  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya completa)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed}" | tee -a "$LOG"

    ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
      sa irandomk tmaxiter 200 nangshift 5 sa_metropolis 50000 1 250 \
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

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE" | tee -a "$LOG"
