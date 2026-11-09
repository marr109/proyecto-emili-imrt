#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_tabu_organs.sh <step: 5|10> <strategy: first|best> [seed_from=1] [seed_to=15]
# Tabu search (token tabubao: neighbours from the CURRENT solution; the framework 'tabu' token only
# explores around the best solution and stalls) with adaptive memory (tenure 3..8), restricted catalog, K=4, per-organ FMO model.
# Terminates by wall-clock budget (BUDGET_SECONDS, default 1080 s = 18 min) via "-it".
# Output: ${OUT_ROOT:-experiments/new}/tabu/nangshift<step>/restricted/data/<strategy>/seedNN/

STEP="$1"; STRATEGY="$2"; SEED_FROM="${3:-1}"; SEED_TO="${4:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

BIN="${EMILI_BIN:-./build/emili}"
BUDGET="${BUDGET_SECONDS:-1080}"
PY="${PY:-$HOME/.venvs/emili-ampl/bin/python}"
BASE="${OUT_ROOT:-experiments/new}/tabu/nangshift${STEP}/restricted"
OUTBASE="${BASE}/data/${STRATEGY}"
LOG="${BASE}/${STRATEGY}_batch.log"
mkdir -p "$OUTBASE"

for n in $(seq "$SEED_FROM" "$SEED_TO"); do
  seed=$(printf '%02d' "$n")
  seed_dir="${OUTBASE}/seed${seed}"; traj="${seed_dir}/trajectory.csv"; run_log="${seed_dir}/run.log"
  mkdir -p "$seed_dir"
  if [ -f "$run_log" ] && grep -q "Conformity Index" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" | tee -a "$LOG"; continue
  fi
  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (tabu ${STRATEGY} step=${STEP})" | tee -a "$LOG"
  # -it must come BEFORE rnds (after rnds it is ignored silently)
  "$BIN" instances/CERR_Prostate baoimrt 4 csv "$traj" \
    tabubao "$STRATEGY" irandomk tmaxiter 1000000 nangshift "$STEP" TBao_adaptive 3 8 \
    -it "$BUDGET" rnds "$n" > "$run_log" 2>&1 || true
  if grep -q "Conformity Index" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no 'Conformity Index', check run.log" | tee -a "$LOG"; continue
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (tabu ${STEP} ${STRATEGY})" | tee -a "$LOG"
