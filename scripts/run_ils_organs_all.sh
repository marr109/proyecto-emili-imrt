#!/usr/bin/env bash
set -euo pipefail

# Runs the whole ILS matrix with the per-organ model, every run with the 1080 s budget:
#   5 perturbations x step{5,10} x catalog{restricted,unrestricted} x
#   { first x order{fixed,circular,rand} + best x fixed } x 15 seeds
#   = 5 x 2 x 2 x 4 = 80 conditions, 1200 runs (best ignores the scan order, so only fixed is run).
# Usage: scripts/run_ils_organs_all.sh [parallel_jobs] [split=1]
# split = workers per condition; the 15 seeds are divided into that many ranges.
# All runs last the same wall-clock budget, so jobs are queued in a fixed order (first, then best).
# Each AMPL/Gurobi session needs a license seat: check how many concurrent sessions your license allows.
# The launcher must export EMILI_GUROBI_OPTIONS="threads=1".
JOBS="${1:-1}"; SPLIT="${2:-1}"
for pert in range3 adaptive mag2 mag2-stepwise adaptivemag-stepwise; do
for step in 5 10; do for catalog in restricted unrestricted; do
for so in "first fixed" "first circular" "first rand" "best fixed"; do
  set -- $so; strategy="$1"; order="$2"
  for part in $(seq 0 $((SPLIT - 1))); do
    echo "$pert $order $step $strategy $catalog $((part * 15 / SPLIT + 1)) $(((part + 1) * 15 / SPLIT))"
  done
done; done; done; done | xargs -P "$JOBS" -L 1 scripts/run_ils_organs.sh
