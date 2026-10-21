#!/usr/bin/env bash
set -euo pipefail

# Runs the whole VND matrix (3 configurations x 15 seeds = 45 runs), one job per (configuration, seed).
# Usage: scripts/run_vnd_organs_all.sh [parallel_jobs]
# Each AMPL/Gurobi session needs a license seat. Keep jobs <= 48 (physical cores of the server).
JOBS="${1:-1}"
for seed in $(seq 1 15); do
  for config in 10-swap 10-5-swap 10-5; do      # slowest configurations first
    echo "$config $seed $seed"
  done
done | xargs -P "$JOBS" -L 1 scripts/run_vnd_organs.sh
