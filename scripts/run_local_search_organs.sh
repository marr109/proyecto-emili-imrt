#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_local_search_organs.sh <step> <strategy: first|best> \
#            <catalog: restricted|unrestricted> <order: fixed|circular|rand> [seed_from=1] [seed_to=15]
# Pure local search (no ILS), 15 seeds, K=4, with the PER-ORGAN FMO model.
# Output (same layout as experiments/old/local_search):
#   experiments/new/local_search/nangshift<step>/<catalog>[-shift-circular|-shift-rand]/data/<strategy>/seedNN/
# Never writes to experiments/old.
#
# Model settings (override via environment):
#   EMILI_FMO_MODEL=organs  EMILI_W_OVER_ORGANS=0.5,0.8 (bladder,rectum)
#   EMILI_DMIN=68  EMILI_W_PTV_OVER=0.5  EMILI_EXCLUDE_OVERLAP=1  EMILI_MAX_INTENSITY=100

STEP="$1"; STRATEGY="$2"; CATALOG="$3"; ORDER="$4"
SEED_FROM="${5:-1}"; SEED_TO="${6:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

case "$CATALOG" in
  unrestricted) CATALOG_TOKEN="unrestricted" ;;
  restricted)   CATALOG_TOKEN="" ;;
  *) echo "catalog must be restricted|unrestricted, got: $CATALOG" >&2; exit 1 ;;
esac
case "$ORDER" in
  fixed)    ORDER_TOKEN="";          SUFFIX="" ;;
  circular) ORDER_TOKEN="circular";  SUFFIX="-shift-circular" ;;
  rand)     ORDER_TOKEN="randorder"; SUFFIX="-shift-rand" ;;
  *) echo "order must be fixed|circular|rand, got: $ORDER" >&2; exit 1 ;;
esac

BASE="experiments/new/local_search/nangshift${STEP}/${CATALOG}${SUFFIX}"
OUTBASE="${BASE}/data/${STRATEGY}"
LOG="${BASE}/${STRATEGY}_batch.log"
mkdir -p "$OUTBASE"

for n in $(seq "$SEED_FROM" "$SEED_TO"); do
  seed=$(printf '%02d' "$n")
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"
  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" | tee -a "$LOG"
    continue
  fi

  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} ${STRATEGY} ${CATALOG} ${ORDER})" | tee -a "$LOG"
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" $ORDER_TOKEN \
    rnds "$n" > "$run_log" 2>&1 || true

  if grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no 'Found solution', check run.log" | tee -a "$LOG"
    continue
  fi

  # Same per-seed plots as the old batches (convergence and iterations).
  if [ ! -f "${seed_dir}/convergence.png" ]; then
    "${PY:-ampl_gurobi/.venv/bin/python3}" scripts/plot_trajectory.py \
      --csv "$traj" --skip-trajectory \
      --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
      --convergence-output "${seed_dir}/convergence.png" \
      --iterations-output "${seed_dir}/iterations.png" >> "$run_log" 2>&1 || true
    echo "[$(date '+%H:%M:%S')] PLOT seed${seed}" | tee -a "$LOG"
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} ${STRATEGY} ${CATALOG} ${ORDER})" | tee -a "$LOG"
