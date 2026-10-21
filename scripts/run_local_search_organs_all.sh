#!/usr/bin/env bash
set -euo pipefail

# Runs the 24 local-search conditions (step x catalog x order x strategy) with the
# per-organ model. Usage: scripts/run_local_search_organs_all.sh [parallel_jobs] [split=1]
# split = workers per condition; the 15 seeds are divided into that many ranges.
# Each job is one condition (or one seed range of it), seeds sequential. Each AMPL/Gurobi session needs a
# license seat: check how many concurrent sessions your license allows.
JOBS="${1:-1}"; SPLIT="${2:-1}"
for step in 5 10; do for catalog in restricted unrestricted; do
for order in fixed circular rand; do for strategy in first best; do
  for part in $(seq 0 $((SPLIT - 1))); do
    echo "$step $strategy $catalog $order $((part * 15 / SPLIT + 1)) $(((part + 1) * 15 / SPLIT))"
  done
done; done; done; done | xargs -P "$JOBS" -L 1 scripts/run_local_search_organs.sh
