#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_sa_geometric_organs.sh <step: 5|10> \
#            <catalog: restricted|unrestricted> [seed_from=1] [seed_to=15]
# Simulated annealing (Metropolis) with the PER-ORGAN FMO model, 15 seeds, K=4.
# Output: ${OUT_ROOT:-experiments/new}/sa/nangshift<step>/<catalog>-budget-geometric/data/seedNN/
# Never writes to experiments/old.
#
# GEOMETRIC-cooling twin of run_sa_organs.sh (same budget, neighborhood, T range, model).
# Uses the framework's 'saacc <start> <end> <beta> <interval> <alpha>' acceptance, where the
# temperature is updated every <interval> calls as T <- alpha*T - beta (emilibase.cpp, Metropolis::accept)
# and cooling stops once T <= end. With beta=0, interval=1 and alpha=0.977 the temperature goes
# from 1300 to 80 in ln(80/1300)/ln(0.977) ~ 120 iterations (the budget allows ~110-135).
#
# Model settings (override via environment):
#   EMILI_FMO_MODEL=organs  EMILI_W_OVER_ORGANS=0.5,0.8 (bladder,rectum)
#   EMILI_DMIN=68  EMILI_W_PTV_OVER=0.5  EMILI_EXCLUDE_OVERLAP=1  EMILI_MAX_INTENSITY=100
# Run-control overrides (for smoke tests):
#   OUT_ROOT (default experiments/new)  BUDGET_SECONDS (default 1080; smoke tests only)
#   T_START/T_END/ALPHA  PY

STEP="$1"; CATALOG="$2"
SEED_FROM="${3:-1}"; SEED_TO="${4:-15}"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

OUT_ROOT="${OUT_ROOT:-experiments/new}"
BUDGET="${BUDGET_SECONDS:-1080}"
T_START="${T_START:-1300}"; T_END="${T_END:-80}"; ALPHA="${ALPHA:-0.977}"

case "$STEP" in 5|10) ;; *) echo "step must be 5|10, got: $STEP" >&2; exit 1 ;; esac
case "$CATALOG" in
  unrestricted) CATALOG_TOKEN="unrestricted" ;;
  restricted)   CATALOG_TOKEN="" ;;
  *) echo "catalog must be restricted|unrestricted, got: $CATALOG" >&2; exit 1 ;;
esac
# -it must come BEFORE rnds (rnds uses seek, -it uses checkToken: after rnds it is ignored silently).
# Budget cuts skip "Found solution"; the clinical report (atexit) is the completion mark.
DONE_MARK="Conformity Index"
RUN_TOKENS="tmaxiter 1000000 nangshift ${STEP} saacc ${T_START} ${T_END} 0 1 ${ALPHA} -it ${BUDGET}"

BASE="${OUT_ROOT}/sa/nangshift${STEP}/${CATALOG}-budget-geometric"
OUTBASE="${BASE}/data"
LOG="${BASE}/batch.log"
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

  echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} ${CATALOG})" | tee -a "$LOG"
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
    sa irandomk $CATALOG_TOKEN $RUN_TOKENS \
    rnds "$n" > "$run_log" 2>&1 || true

  if grep -q "$DONE_MARK" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: no '${DONE_MARK}', check run.log" | tee -a "$LOG"
    continue
  fi

  if [ ! -f "${seed_dir}/convergence.png" ]; then
    "${PY:-ampl_gurobi/.venv/bin/python3}" scripts/plot_trajectory.py \
      --csv "$traj" --skip-trajectory \
      --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
      --convergence-output "${seed_dir}/convergence.png" \
      --iterations-output "${seed_dir}/iterations.png" >> "$run_log" 2>&1 || true
    echo "[$(date '+%H:%M:%S')] PLOT seed${seed}" | tee -a "$LOG"
  fi
done
echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} ${CATALOG})" | tee -a "$LOG"
