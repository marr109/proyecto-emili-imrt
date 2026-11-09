#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_vns_organs.sh <config> [seed_from=1] [seed_to=15]
# General VNS (shake + neighbourhood change), first improvement, restricted catalog, K=4, per-organ FMO model.
# Configs (all: acceptance 'baoimprove rejectrepeated'):
#   rvns         : local search shift 10, multi-scale BAO shake (levels up to 3)
#   pershake123  : local search shift 10, shake = perturbations of 1, 2 and 3 angles (shift 10)
#   vnd-rvns     : VND (shift 10 then 5) as local search, multi-scale BAO shake
#   vnd-swap     : VND (shift 10, 5, sampled swap of 10 neighbours) as local search, multi-scale BAO shake
# Wall-clock budget BUDGET_SECONDS (default 1500). Analysis truncates every run at a common number of real FMO solves.
# EMILI_BIN overrides the binary (default ./build/emili). Output: ${OUT_ROOT:-experiments/new}/vns/<config>/data/first/seedNN/

CONFIG="$1"; SEED_FROM="${2:-1}"; SEED_TO="${3:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

BIN="${EMILI_BIN:-./build/emili}"
BUDGET="${BUDGET_SECONDS:-1500}"
ACC="accng baoimprove rejectrepeated"
case "$CONFIG" in
  rvns)        ALG="vns first irandomk locmin nangshift 10";                            SHAKE="bangshake 3" ;;
  pershake123) ALG="vns first irandomk locmin nangshift 10";                            SHAKE="pershake prangshift 10 1 prangshift 10 2 prangshift 10 3" ;;
  vnd-rvns)    ALG="vns vnd first irandomk locmin nangshift 10 nangshift 5";            SHAKE="bangshake 3" ;;
  vnd-swap)    ALG="vns vnd first irandomk locmin nangshift 10 nangshift 5 nangswapsample 10"; SHAKE="bangshake 3" ;;
  *) echo "unknown config: $CONFIG" >&2; exit 1 ;;
esac

BASE="${OUT_ROOT:-experiments/new}/vns/${CONFIG}"
OUTBASE="${BASE}/data/first"
LOG="${BASE}/first_batch.log"
mkdir -p "$OUTBASE"

for n in $(seq "$SEED_FROM" "$SEED_TO"); do
  seed=$(printf '%02d' "$n")
  seed_dir="${OUTBASE}/seed${seed}"; traj="${seed_dir}/trajectory.csv"; run_log="${seed_dir}/run.log"
  mkdir -p "$seed_dir"
  if [ -f "$run_log" ] && grep -q "Conformity Index" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" | tee -a "$LOG"; continue
  fi
  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (${CONFIG})" | tee -a "$LOG"
  # -it must come BEFORE rnds (after rnds it is ignored silently)
  "$BIN" instances/CERR_Prostate baoimrt 4 csv "$traj" \
    $ALG tmaxiter 1000000 $SHAKE $ACC -it "$BUDGET" rnds "$n" > "$run_log" 2>&1 || true
  if grep -q "Conformity Index" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no 'Conformity Index', check run.log" | tee -a "$LOG"
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (${CONFIG})" | tee -a "$LOG"
