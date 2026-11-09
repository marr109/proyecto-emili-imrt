#!/usr/bin/env bash
set -euo pipefail
# Usage: scripts/run_vns_organs_all.sh <parallel_jobs> <config> [config...]   (one job per config and seed)
JOBS="$1"; shift
for seed in $(seq 1 15); do for cfg in "$@"; do echo "$cfg $seed $seed"; done; done | xargs -P "$JOBS" -L 1 scripts/run_vns_organs.sh
