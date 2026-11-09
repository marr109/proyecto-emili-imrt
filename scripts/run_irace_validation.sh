#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_irace_validation.sh <config> <seed>
# Validation of the elite configurations found by irace (irace/results/irace.Rdata, final elites 43, 45, 33).
# One run = one (config, seed) pair, per-organ FMO model, K=4, wall-clock budget BUDGET_SECONDS (default 1080 = 18 min).
# Seeds 1..15 are the ones used by the baselines and differ from the irace seeds.
# Configs (all: VNS, VND as local search, first improvement, acceptance 'accng baoimprove rejectrepeated'):
#   irace43 : sampled swap (18) + window swap (10 deg); shake = magnitude perturbation of 2 angles, stepwise
#   irace45 : sampled swap (24) + window swap (45 deg); shake = magnitude perturbation of 2 angles
#   irace33 : sampled swap (16) + window swap (25 deg) + flip; shake = multi-scale BAO shake (levels up to 4)
# Output: ${OUT_ROOT:-experiments/new/irace_validation}/<config>/data/first/seedNN/

CONFIG="$1"; N="$2"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

BIN="${EMILI_BIN:-./build/emili}"
BUDGET="${BUDGET_SECONDS:-1080}"
ACC="accng baoimprove rejectrepeated"
case "$CONFIG" in
  irace43) ALG="vns vnd first irandomk locmin nangswapsample 18 nangswapwindow 10";               SHAKE="pershake prangshiftmag 2 2 5 10 stepwise" ;;
  irace45) ALG="vns vnd first irandomk locmin nangswapsample 24 nangswapwindow 45";               SHAKE="pershake prangshiftmag 2 2 5 10" ;;
  irace33) ALG="vns vnd first irandomk locmin nangswapsample 16 nangswapwindow 25 nangflip 0";    SHAKE="bangshake 4" ;;
  # Hand-designed reference (expanded VND) under the same 18-min budget and machine load as the elites
  vndref)  ALG="vnd first irandomk locmin nangshift 10 nangshift 5 nangswapsample 10 nangswapwindow 20 nangflip 0"; SHAKE=""; ACC="" ;;
  *) echo "unknown config: $CONFIG" >&2; exit 1 ;;
esac
TAIL_ARGS="tmaxiter 1000000 $SHAKE $ACC"
[ "$CONFIG" = "vndref" ] && TAIL_ARGS=""

BASE="${OUT_ROOT:-experiments/new/irace_validation}/${CONFIG}"
seed=$(printf '%02d' "$N")
seed_dir="${BASE}/data/first/seed${seed}"; traj="${seed_dir}/trajectory.csv"; run_log="${seed_dir}/run.log"
LOG="${BASE}/first_batch.log"
mkdir -p "$seed_dir"
if [ -f "$run_log" ] && grep -q "Conformity Index" "$run_log" 2>/dev/null; then
  echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" >> "$LOG"; exit 0
fi
echo "[$(date '+%H:%M:%S')] RUN seed${seed} (${CONFIG})" >> "$LOG"
# -it must come BEFORE rnds (after rnds it is ignored silently)
"$BIN" instances/CERR_Prostate baoimrt 4 csv "$traj" \
  $ALG $TAIL_ARGS -it "$BUDGET" rnds "$N" > "$run_log" 2>&1 || true
if grep -q "Conformity Index" "$run_log" 2>/dev/null; then
  echo "[$(date '+%H:%M:%S')] DONE seed${seed} (${CONFIG})" >> "$LOG"
else
  echo "[$(date '+%H:%M:%S')] WARNING seed${seed} (${CONFIG}): no 'Conformity Index', check run.log" >> "$LOG"
fi
