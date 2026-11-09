#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/run_catalog_job.sh <family> <config> <catalog: restricted|unrestricted> <seed>
# One run = one (family, config, catalog, seed), per-organ FMO model, K=4, wall-clock budget BUDGET_SECONDS
# (default 1080 = 18 min). "unrestricted" adds the token to irandomk and to every swap/window/flip
# neighbourhood (nangshift has no catalog token: it moves +-step degrees from the current angle).
#
# Families and configs:
#   vnd      10-5 | 10-5-sample | 10-5-sample-window-flip   (VND, first improvement; ends at the local optimum,
#                                                            safety limit SAFETY_SECONDS, default 5400)
#   vns      rvns | pershake123 | vnd-rvns | vnd-swap       (general VNS, as scripts/run_vns_organs.sh)
#   tabu     <step 5|10>-<first|best>                       (tabubao, TBao_adaptive 3 8, as run_tabu_organs.sh)
#   tabuattr t5 | t10 | t20 | t10-5-first | t10-10-best | t10-5-best   (tabubao <strategy>, nangshift <step>, TBao_attr <tenure>)
#
# Output (the restricted layouts of the existing batches are not touched):
#   vnd      ${OUT_ROOT:-experiments/new}/vnd/nangshift<config>-unrestricted/data/first/seedNN/
#   vns      ${OUT_ROOT}/vns/<config>-unrestricted/data/first/seedNN/
#   tabu     ${OUT_ROOT}/tabu/nangshift<step>/unrestricted/data/<strategy>/seedNN/
#   tabuattr ${OUT_ROOT}/tabu_attr/nangshift10-<config>/<catalog>/data/first/seedNN/

FAMILY="$1"; CONFIG="$2"; CATALOG="$3"; N="$4"

export EMILI_FMO_MODEL="${EMILI_FMO_MODEL:-organs}"
export EMILI_W_OVER_ORGANS="${EMILI_W_OVER_ORGANS:-0.5,0.8}"
export EMILI_DMIN="${EMILI_DMIN:-68}"
export EMILI_W_PTV_OVER="${EMILI_W_PTV_OVER:-0.5}"
export EMILI_EXCLUDE_OVERLAP="${EMILI_EXCLUDE_OVERLAP:-1}"
export EMILI_MAX_INTENSITY="${EMILI_MAX_INTENSITY:-100}"

# Linux: resolve the AMPL / Gurobi binaries of the repo venv once (see run_ils_budget_batch.sh)
if [ -z "${EMILI_AMPL_BIN_DIR:-}" ] || [ -z "${EMILI_GUROBI_BIN:-}" ]; then
  SP="$(echo ampl_gurobi/.venv/lib/python*/site-packages)"
  if [ -d "$SP/ampl_module_base" ]; then
    export EMILI_AMPL_BIN_DIR="$PWD/$SP/ampl_module_base/bin"
    export EMILI_GUROBI_BIN="$PWD/$SP/ampl_module_gurobi/bin/gurobi"
  fi
fi

BIN="${EMILI_BIN:-./build/emili}"
BUDGET="${BUDGET_SECONDS:-1080}"
SAFETY="${SAFETY_SECONDS:-5400}"

case "$CATALOG" in
  unrestricted) U=" unrestricted" ;;
  restricted)   U="" ;;
  *) echo "catalog must be restricted|unrestricted, got: $CATALOG" >&2; exit 1 ;;
esac
SAMPLE="nangswapsample 10${U}"; WINDOW="nangswapwindow 20${U}"; FLIP="nangflip 0${U}"
TIME_LIMIT="$BUDGET"

case "$FAMILY" in
  vnd)
    case "$CONFIG" in
      10-5)                    NEIGH="nangshift 10 nangshift 5" ;;
      10-5-sample)             NEIGH="nangshift 10 nangshift 5 $SAMPLE" ;;
      10-5-sample-window-flip) NEIGH="nangshift 10 nangshift 5 $SAMPLE $WINDOW $FLIP" ;;
      *) echo "unknown vnd config: $CONFIG" >&2; exit 1 ;;
    esac
    ARGS="vnd first irandomk${U} locmin $NEIGH"; TIME_LIMIT="$SAFETY"
    BASE="${OUT_ROOT:-experiments/new}/vnd/nangshift${CONFIG}-${CATALOG}"; SUB="first" ;;
  vns)
    ACC="accng baoimprove rejectrepeated"
    case "$CONFIG" in
      rvns)        ALG="vns first irandomk${U} locmin nangshift 10";                         SHAKE="bangshake 3" ;;
      pershake123) ALG="vns first irandomk${U} locmin nangshift 10";                         SHAKE="pershake prangshift 10 1 prangshift 10 2 prangshift 10 3" ;;
      vnd-rvns)    ALG="vns vnd first irandomk${U} locmin nangshift 10 nangshift 5";         SHAKE="bangshake 3" ;;
      vnd-swap)    ALG="vns vnd first irandomk${U} locmin nangshift 10 nangshift 5 $SAMPLE"; SHAKE="bangshake 3" ;;
      *) echo "unknown vns config: $CONFIG" >&2; exit 1 ;;
    esac
    ARGS="$ALG tmaxiter 1000000 $SHAKE $ACC"
    BASE="${OUT_ROOT:-experiments/new}/vns/${CONFIG}-${CATALOG}"; SUB="first" ;;
  tabu)
    STEP="${CONFIG%-*}"; STRATEGY="${CONFIG#*-}"
    ARGS="tabubao $STRATEGY irandomk${U} tmaxiter 1000000 nangshift $STEP TBao_adaptive 3 8"
    BASE="${OUT_ROOT:-experiments/new}/tabu/nangshift${STEP}/${CATALOG}"; SUB="$STRATEGY" ;;
  tabuattr)
    # config: t<tenure>[-<step>[-<strategy>]]  (defaults: step 10, first)
    IFS=- read -r TENURE STEP STRATEGY <<< "${CONFIG#t}"
    STEP="${STEP:-10}"; STRATEGY="${STRATEGY:-first}"
    case "$STEP" in 1|5|10|20) ;; *) echo "tabuattr: bad step '$STEP' in config '$CONFIG'" >&2; exit 1 ;; esac
    case "$STRATEGY" in first|best) ;; *) echo "tabuattr: bad strategy '$STRATEGY' in config '$CONFIG'" >&2; exit 1 ;; esac
    ARGS="tabubao $STRATEGY irandomk${U} tmaxiter 1000000 nangshift $STEP TBao_attr $TENURE"
    BASE="${OUT_ROOT:-experiments/new}/tabu_attr/nangshift${STEP}-t${TENURE}/${CATALOG}"; SUB="$STRATEGY" ;;
  *) echo "unknown family: $FAMILY" >&2; exit 1 ;;
esac

seed=$(printf '%02d' "$N")
seed_dir="${BASE}/data/${SUB}/seed${seed}"; traj="${seed_dir}/trajectory.csv"; run_log="${seed_dir}/run.log"
LOG="${BASE}/batch.log"
mkdir -p "$seed_dir"
if [ -f "$run_log" ] && grep -q "Conformity Index" "$run_log" 2>/dev/null; then
  echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (already complete)" >> "$LOG"; exit 0
fi
echo "[$(date '+%H:%M:%S')] RUN seed${seed} (${FAMILY} ${CONFIG} ${CATALOG})" >> "$LOG"
t0=$(date +%s)
# -it must come BEFORE rnds (after rnds it is ignored silently)
"$BIN" "${INSTANCE:-instances/CERR_Prostate}" baoimrt 4 csv "$traj" \
  $ARGS -it "$TIME_LIMIT" rnds "$N" > "$run_log" 2>&1 || true
elapsed=$(( $(date +%s) - t0 ))
if grep -q "Conformity Index" "$run_log" 2>/dev/null; then
  note=""
  [ "$FAMILY" = "vnd" ] && [ "$elapsed" -ge $(( SAFETY - 10 )) ] && note=" HIT-SAFETY-LIMIT"
  echo "[$(date '+%H:%M:%S')] DONE seed${seed} (${FAMILY} ${CONFIG} ${CATALOG}) ${elapsed}s${note}" >> "$LOG"
else
  echo "[$(date '+%H:%M:%S')] WARNING seed${seed} (${FAMILY} ${CONFIG} ${CATALOG}): no 'Conformity Index', check run.log" >> "$LOG"
fi
