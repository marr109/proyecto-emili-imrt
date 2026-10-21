#!/usr/bin/env bash
# =============================================================================
# Comparación final de metaheurísticas BAO en PROSTATE_36ang K=4
# Cada algoritmo corre con sus parámetros propios en paralelo.
# OSQP ~240s/call en esta instancia → runs de 40-90 min esperados.
# IG excluido: reconstrucción greedy = 33 × 240s = 2.2h por perturbación.
# =============================================================================
set -euo pipefail

EMILI="./build/emili"
INST="instances/PROSTATE_36ang"
K=4
SEED=42
OUTDIR="comparison_results"
mkdir -p "$OUTDIR"

run_algo() {
    local name="$1"; shift
    local logfile="$OUTDIR/${name}.log"
    echo "[$(date '+%H:%M:%S')] Iniciando $name..."
    time "$EMILI" "$INST" baoimrt "$K" "$@" rnds "$SEED" > "$logfile" 2>&1
    f=$(grep "Objective function value:" "$logfile" | awk '{print $NF}')
    angles=$(grep "Found solution:" "$logfile" | head -1)
    echo "[$(date '+%H:%M:%S')] $name → f=$f"
    echo "  $angles"
}

# Lanzar en paralelo (4 algoritmos, 1 core cada uno)
run_algo "ILS" \
    ils first irandomk tmaxiter 1 nangswap tmaxiter 3 prangswap 2 improve &

run_algo "VNS" \
    vns first irandomk tmaxiter 1 nangswap tmaxiter 2 bangshake 2 accng improve &

run_algo "SA" \
    ils nols irandomk tmaxiter 15 prangswap 1 \
    saacc 5000 50 0.1 1 0.95 &

run_algo "Tabu" \
    tabu first irandomk tmaxiter 5 nangswap TBao_fixed 3 &

echo "4 algoritmos corriendo en paralelo. Esperando resultados..."
wait

echo ""
echo "===== RESULTADOS FINALES ====="
for algo in ILS VNS SA Tabu; do
    f=$(grep "Objective function value:" "$OUTDIR/${algo}.log" | awk '{print $NF}')
    angles=$(grep "Found solution:" "$OUTDIR/${algo}.log" | head -1 | sed 's/Found solution: //')
    echo "$algo:  f=$f  |  $angles"
done
