#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/launch_catalog_study.sh [jobs_file]   (run from the repo root, detached with nohup)
# Catalog study (restricted vs unrestricted) + attribute tabu memory: builds the job queue in priority order
# and runs it with at most PARALLEL (default 48) concurrent runs, one Gurobi thread each, low priority.
#   1. tabu with attribute memory (TBao_attr 5/10/20), restricted and unrestricted      90 runs
#   2. general VNS unrestricted (rvns, pershake123, vnd-rvns, vnd-swap)                 60 runs
#   3. tabu with TBao_adaptive unrestricted (step 5/10 x first/best)                    60 runs
#   4. VND unrestricted (10-5, 10-5-sample, 10-5-sample-window-flip)                    45 runs
# Each job is "family config catalog seed"; scripts/run_catalog_job.sh skips runs that are already complete,
# so the launcher can be re-run to resume.

JOBS="${1:-experiments/new/catalog_study_jobs.txt}"
PARALLEL="${PARALLEL:-48}"
export EMILI_GUROBI_OPTIONS="${EMILI_GUROBI_OPTIONS:-threads=1}"
mkdir -p "$(dirname "$JOBS")"
: > "$JOBS"
add() { local fam="$1" catalog="$2"; shift 2
  for cfg in "$@"; do for s in $(seq 1 15); do echo "$fam $cfg $catalog $s" >> "$JOBS"; done; done; }
add tabuattr restricted   t5 t10 t20
add tabuattr unrestricted t5 t10 t20
add vns      unrestricted rvns pershake123 vnd-rvns vnd-swap
add tabu     unrestricted 10-first 10-best 5-first 5-best
add vnd      unrestricted 10-5 10-5-sample 10-5-sample-window-flip
echo "[$(date '+%F %T')] queue: $(wc -l < "$JOBS") jobs, parallel=$PARALLEL, $EMILI_GUROBI_OPTIONS"
xargs -P "$PARALLEL" -L 1 nice -n 10 scripts/run_catalog_job.sh < "$JOBS"
echo "[$(date '+%F %T')] QUEUE COMPLETE"
