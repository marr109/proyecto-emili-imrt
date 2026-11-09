#!/usr/bin/env bash
set -euo pipefail

# Usage: OUT_ROOT=<dir> scripts/launch_catalog_restricted.sh [jobs_file]   (run from the repo root, detached with nohup)
# Restricted-catalog counterpart of scripts/launch_catalog_study.sh, run on the same machine as the unrestricted batch so the
# two catalogs can be compared without a platform effect. PARALLEL defaults to 32: on the 48-core node, 48 concurrent runs do
# ~72 FMO solves per 18 minutes, 32 do ~90 (the Mac does ~92). Each job is "family config catalog seed".
#   1. general VNS restricted (rvns, pershake123, vnd-rvns, vnd-swap)                 60 runs
#   2. tabu with TBao_adaptive restricted (step 5/10 x first/best)                     60 runs
#   3. VND restricted (10-5, 10-5-sample, 10-5-sample-window-flip)                     45 runs

JOBS="${1:-experiments/new/catalog_restricted_jobs.txt}"
PARALLEL="${PARALLEL:-32}"
export EMILI_GUROBI_OPTIONS="${EMILI_GUROBI_OPTIONS:-threads=1}"
mkdir -p "$(dirname "$JOBS")"
: > "$JOBS"
add() { local fam="$1" catalog="$2"; shift 2
  for cfg in "$@"; do for s in $(seq 1 15); do echo "$fam $cfg $catalog $s" >> "$JOBS"; done; done; }
add vns  restricted rvns pershake123 vnd-rvns vnd-swap
add tabu restricted 10-first 10-best 5-first 5-best
add vnd  restricted 10-5 10-5-sample 10-5-sample-window-flip
echo "[$(date '+%F %T')] queue: $(wc -l < "$JOBS") jobs, parallel=$PARALLEL, $EMILI_GUROBI_OPTIONS"
xargs -P "$PARALLEL" -L 1 nice -n 10 scripts/run_catalog_job.sh < "$JOBS"
echo "[$(date '+%F %T')] QUEUE COMPLETE"
