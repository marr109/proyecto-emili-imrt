#!/usr/bin/env bash
set -euo pipefail

# Runs the whole GEOMETRIC-cooling SA matrix with the per-organ model: step{5,10} x catalog
# = 4 conditions x 15 seeds = 60 runs.
# Usage: scripts/run_sa_organs_all.sh [parallel_jobs] [split=1]
# split = workers per condition; the 15 seeds are divided into that many ranges.
# The launcher must export EMILI_GUROBI_OPTIONS="threads=1".
JOBS="${1:-1}"; SPLIT="${2:-1}"
for step in 5 10; do for catalog in restricted unrestricted; do
  for part in $(seq 0 $((SPLIT - 1))); do
    echo "$step $catalog $((part * 15 / SPLIT + 1)) $(((part + 1) * 15 / SPLIT))"
  done
done; done | xargs -P "$JOBS" -L 1 scripts/run_sa_geometric_organs.sh
