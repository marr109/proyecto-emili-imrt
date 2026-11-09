#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/regen_seed_plots.sh <experiments_subdir>
# Example: scripts/regen_seed_plots.sh experiments/new/ils
# Regenerates convergence.png, iterations.png and trajectory_improvements.csv
# in every seedNN directory (containing trajectory.csv) under <experiments_subdir>.
# These files are not kept in the repository layout (only run.log and trajectory.csv).

if [ $# -ne 1 ] || [ ! -d "$1" ]; then
  echo "Usage: $0 <experiments_subdir>" >&2
  exit 1
fi

PY="${PY:-$HOME/.venvs/emili-ampl/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

count=0
while IFS= read -r traj; do
  seed_dir="$(dirname "$traj")"
  "$PY" "${SCRIPT_DIR}/plot_trajectory.py" \
    --csv "$traj" \
    --skip-trajectory \
    --improvements-csv "${seed_dir}/trajectory_improvements.csv" \
    --convergence-output "${seed_dir}/convergence.png" \
    --iterations-output "${seed_dir}/iterations.png" > /dev/null
  count=$((count + 1))
done < <(fd --no-ignore -t f -p '/seed[0-9]{2}/trajectory\.csv$' "$1" | sort)
echo "Regenerated plots for ${count} seed directories under $1"
