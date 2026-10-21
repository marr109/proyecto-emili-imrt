#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_vnd_organs.sh <config: 10-5|10-swap|10-5-swap> [seed_from=1] [seed_to=15]
# VND (variable neighborhood descent, first improvement, fixed scan order) with the PER-ORGAN FMO model.
# Restricted catalog, K=4, 15 seeds; terminates by convergence (locmin), no time budget.
#   10-5       : nangshift 10, then nangshift 5
#   10-swap    : nangshift 10, then nangswap (replace one angle by an inactive one)
#   10-5-swap  : nangshift 10, nangshift 5, then nangswap
# Output: ${OUT_ROOT:-experiments/new}/vns/nangshift<config>/data/first/seedNN/{trajectory.csv,run.log}
#
# Model settings (override via environment): same defaults as the other per-organ batches.

CONFIG="$1"; SEED_FROM="${2:-1}"; SEED_TO="${3:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

case "$CONFIG" in
  10-5)      NEIGH="nangshift 10 nangshift 5" ;;
  10-swap)   NEIGH="nangshift 10 nangswap" ;;
  10-5-swap) NEIGH="nangshift 10 nangshift 5 nangswap" ;;
  *) echo "config must be 10-5|10-swap|10-5-swap, got: $CONFIG" >&2; exit 1 ;;
esac

BASE="${OUT_ROOT:-experiments/new}/vns/nangshift${CONFIG}"
OUTBASE="${BASE}/data/first"
LOG="${BASE}/first_batch.log"
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

  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (vnd first ${NEIGH})" | tee -a "$LOG"
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    vnd first irandomk locmin $NEIGH \
    rnds "$n" > "$run_log" 2>&1 || true

  if grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no 'Found solution', check run.log" | tee -a "$LOG"
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (vnd ${CONFIG})" | tee -a "$LOG"
