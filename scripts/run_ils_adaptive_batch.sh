#!/usr/bin/env bash
set -euo pipefail

# Uso: scripts/run_ils_adaptive_batch.sh <step> <first|best> <restricted|unrestricted> [maxNumSteps] [stagnationThreshold]
# Igual que run_ils_batch.sh, pero con perturbación adaptativa progresiva
# (AdaptiveAngleShiftPerturbation, token prangshiftadaptive <step>
# <maxNumSteps> <stagnationThreshold>) en vez de la perturbación fija de 3
# ángulos. Arranca moviendo 1 ángulo (nivel 0) y escala de a uno hasta
# maxNumSteps tras stagnationThreshold rondas seguidas sin mejora del
# objetivo; vuelve a nivel 0 con una mejora real. La magnitud del
# desplazamiento queda fija en (step, 2*step) en todos los niveles.
#
# maxNumSteps y stagnationThreshold son opcionales, default 4 y 2
# (verificados con smoke test real). Salida en:
#   experiments/ils/nangshift<step>/<catalog>-adaptive/data/<strategy>/seedNN/
#
# NO toca las carpetas de orden fijo (grid5/, restricted/, unrestricted/)
# ni las de circular.

STEP="$1"
STRATEGY="$2"
CATALOG="$3"
MAXNUMSTEPS="${4:-4}"
STAGNATION="${5:-2}"
PY="ampl_gurobi/.venv/bin/python3"

if [ "$CATALOG" = "unrestricted" ]; then
  CATALOG_TOKEN="unrestricted"
elif [ "$CATALOG" = "restricted" ]; then
  CATALOG_TOKEN=""
else
  echo "catalog debe ser 'restricted' o 'unrestricted', recibido: $CATALOG" >&2
  exit 1
fi

OUTBASE="experiments/ils/nangshift${STEP}/${CATALOG}-adaptive/data/${STRATEGY}"
LOG="experiments/ils/nangshift${STEP}/${CATALOG}-adaptive/${STRATEGY}_batch.log"

mkdir -p "$OUTBASE"

for seed in $(seq -w 1 15); do
  seed_dir="${OUTBASE}/seed${seed}"
  traj="${seed_dir}/trajectory.csv"
  run_log="${seed_dir}/run.log"

  mkdir -p "$seed_dir"

  if [ -f "$run_log" ] && grep -q "Found solution" "$run_log" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP seed${seed} (ya completa)" | tee -a "$LOG"
  else
    echo "[$(date '+%H:%M:%S')] RUN seed${seed} (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} adaptive)" | tee -a "$LOG"

    ./build/emili instances/CERR_Prostate baoimrt 4 csv "$traj" \
      ils "$STRATEGY" irandomk $CATALOG_TOKEN locmin nangshift "$STEP" \
      tmaxiter 10 \
      prangshiftadaptive "$STEP" "$MAXNUMSTEPS" "$STAGNATION" \
      baoimprove rejectrepeated \
      rnds "$seed" \
      > "$run_log" 2>&1

    if grep -q "Found solution" "$run_log" 2>/dev/null; then
      echo "[$(date '+%H:%M:%S')] DONE seed${seed}" | tee -a "$LOG"
    else
      echo "[$(date '+%H:%M:%S')] WARNING seed${seed}: sin 'Found solution', revisar" | tee -a "$LOG"
      continue
    fi
  fi

  if [ -f "${seed_dir}/convergence.png" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP plot seed${seed}" | tee -a "$LOG"
    continue
  fi

  "$PY" scripts/plot_trajectory.py \
    --csv "$traj" \
    --skip-trajectory \
    --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
    --convergence-output "${seed_dir}/convergence.png" \
    --iterations-output "${seed_dir}/iterations.png" \
    >> "$run_log" 2>&1

  echo "[$(date '+%H:%M:%S')] PLOT seed${seed}" | tee -a "$LOG"
done

echo "[$(date '+%H:%M:%S')] BATCH COMPLETE (step=${STEP} strategy=${STRATEGY} catalog=${CATALOG} adaptive)" | tee -a "$LOG"
