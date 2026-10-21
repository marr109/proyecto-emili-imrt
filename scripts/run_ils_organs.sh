#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_ils_organs.sh <perturbation> <order: fixed|circular|rand> <step: 5|10> \
#            <strategy: first|best> <catalog: restricted|unrestricted> [seed_from=1] [seed_to=15]
# ILS (15 seeds, K=4) with the PER-ORGAN FMO model.
# Output (same layout as experiments/old/ils):
#   ${OUT_ROOT:-experiments/new}/ils/nangshift<step>/<catalog><perturbation-suffix><order-suffix>/data/<strategy>/seedNN/
# Never writes to experiments/old.
#
# EVERY run terminates by a wall-clock budget (BUDGET_SECONDS, default 1080 s = 18 min) via "-it",
# with an unreachable tmaxiter, exactly as run_ils_budget_batch.sh. Extra FMO evaluations of the
# stepwise variants count inside the same budget.
# Perturbation (suffix) -- all with baoimprove rejectrepeated:
#   range3                [ ]                      prangshift <step> 3   (old base + old '-budget1080s')
#   adaptive              [-adaptive]              prangshiftadaptive <step> 4 2
#   mag2                  [-mag2]                  prangshiftmag 2 2 5 10
#   mag2-stepwise         [-mag2-stepwise]         prangshiftmag 2 2 5 10 stepwise
#   adaptivemag-stepwise  [-adaptivemag-stepwise]  prangshiftadaptivemag 2 2 2 5 10 stepwise
# Rule: no new-perturbation variant moves more than 2 angles per perturbation.
# Scan order (suffix): fixed [ ], circular [-shift-circular], rand [-shift-rand] (token randorder).
#   The old '-circular-tmax3' folder maps to range3 + circular (now budget-terminated, not tmaxiter 3).
#   Order only matters for 'first' (best scans the whole neighborhood): best with a non-fixed order is rejected.
#
# Model settings (override via environment):
#   EMILI_FMO_MODEL=organs  EMILI_W_OVER_ORGANS=0.5,0.8 (bladder,rectum)
#   EMILI_DMIN=68  EMILI_W_PTV_OVER=0.5  EMILI_EXCLUDE_OVERLAP=1  EMILI_MAX_INTENSITY=100
# Run-control overrides (for smoke tests):
#   OUT_ROOT (default experiments/new)  BUDGET_SECONDS (default 1080; smoke tests only)  PY

VARIANT="$1"; ORDER="$2"; STEP="$3"; STRATEGY="$4"; CATALOG="$5"
SEED_FROM="${6:-1}"; SEED_TO="${7:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

OUT_ROOT="${OUT_ROOT:-experiments/new}"
BUDGET="${BUDGET_SECONDS:-1080}"

case "$STEP" in 5|10) ;; *) echo "step must be 5|10, got: $STEP" >&2; exit 1 ;; esac
case "$STRATEGY" in first|best) ;; *) echo "strategy must be first|best, got: $STRATEGY" >&2; exit 1 ;; esac
case "$CATALOG" in
  unrestricted) CATALOG_TOKEN="unrestricted" ;;
  restricted)   CATALOG_TOKEN="" ;;
  *) echo "catalog must be restricted|unrestricted, got: $CATALOG" >&2; exit 1 ;;
esac

# -it must come BEFORE rnds (getTime() uses checkToken, rnds uses seek: after rnds -it is ignored silently).
case "$ORDER" in
  fixed)    ORDER_TOKEN="";          OSUF="" ;;
  circular) ORDER_TOKEN="circular";  OSUF="-shift-circular" ;;
  rand)     ORDER_TOKEN="randorder"; OSUF="-shift-rand" ;;
  *) echo "order must be fixed|circular|rand, got: $ORDER" >&2; exit 1 ;;
esac
if [ "$STRATEGY" = best ] && [ "$ORDER" != fixed ]; then
  echo "order has no effect with strategy best; use fixed" >&2; exit 1
fi
case "$VARIANT" in
  range3)               PSUF="";                      PERT="prangshift $STEP 3" ;;
  adaptive)             PSUF="-adaptive";             PERT="prangshiftadaptive $STEP 4 2" ;;
  mag2)                 PSUF="-mag2";                 PERT="prangshiftmag 2 2 5 10" ;;
  mag2-stepwise)        PSUF="-mag2-stepwise";        PERT="prangshiftmag 2 2 5 10 stepwise" ;;
  adaptivemag-stepwise) PSUF="-adaptivemag-stepwise"; PERT="prangshiftadaptivemag 2 2 2 5 10 stepwise" ;;
  *) echo "unknown perturbation: $VARIANT" >&2; exit 1 ;;
esac
SUFFIX="${PSUF}${OSUF}"
# Budget cuts end via finalise()/exit(0), which skips "Found solution"; the clinical report is the completion mark.
DONE_MARK="Conformity Index"

BASE="${OUT_ROOT}/ils/nangshift${STEP}/${CATALOG}${SUFFIX}"
OUTBASE="${BASE}/data/${STRATEGY}"
LOG="${BASE}/${STRATEGY}_batch.log"
mkdir -p "$OUTBASE"

for n in $(seq "$SEED_FROM" "$SEED_TO"); do
  seed=$(printf '%02d' "$n")
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"
  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "$DONE_MARK" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" | tee -a "$LOG"
    continue
  fi

  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (pert=${VARIANT} order=${ORDER} step=${STEP} ${STRATEGY} ${CATALOG})" | tee -a "$LOG"
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    ils "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" $ORDER_TOKEN \
    tmaxiter 1000000 \
    $PERT \
    baoimprove rejectrepeated \
    -it "$BUDGET" \
    rnds "$n" > "$run_log" 2>&1 || true

  if grep -q "$DONE_MARK" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no '${DONE_MARK}', check run.log" | tee -a "$LOG"
    continue
  fi

  # Same per-seed plots as the old batches (convergence and iterations).
  # plot_ils_trajectory.py is NOT called: it re-simulates ILS semantics and is only valid for
  # first + the original perturbations.
  if [ ! -f "${seed_dir}/convergence.png" ]; then
    "${PY:-ampl_gurobi/.venv/bin/python3}" scripts/plot_trajectory.py \
      --csv "$traj" --skip-trajectory \
      --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
      --convergence-output "${seed_dir}/convergence.png" \
      --iterations-output "${seed_dir}/iterations.png" >> "$run_log" 2>&1 || true
    echo "[$(date '+%H:%M:%S')] PLOT seed${seed}" | tee -a "$LOG"
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (pert=${VARIANT} order=${ORDER} step=${STEP} ${STRATEGY} ${CATALOG})" | tee -a "$LOG"
