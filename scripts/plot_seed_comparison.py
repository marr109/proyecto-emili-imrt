#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overlay First vs Best Improvement trajectories for a single seed (same starting
point, same neighborhood): running-best step line over the faint raw evaluation
trace, for each strategy.

CLI example:
  python scripts/plot_seed_comparison.py \
      --first-csv experiments/X/data/first/seed01/trajectory.csv \
      --best-csv experiments/X/data/best/seed01/trajectory.csv \
      --seed 1 --condition-label grid5 \
      --output experiments/X/plots/seed_comparisons/seed01_comparison.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory, running_best


def plot_pair(first_csv, best_csv, seed, condition_label, output):
    f_evals, f_obj, _ = read_trajectory(first_csv)
    b_evals, b_obj, _ = read_trajectory(best_csv)
    f_best = running_best(f_obj)
    b_best = running_best(b_obj)

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(f_evals, f_obj, color="C0", linewidth=0.8, alpha=0.35)
    ax.step(f_evals, f_best, where="post", color="C0", linewidth=2.5,
            label="First (final={0:,.0f}, {1} eval.)".format(f_best[-1], len(f_evals)))

    ax.plot(b_evals, b_obj, color="C1", linewidth=0.8, alpha=0.35)
    ax.step(b_evals, b_best, where="post", color="C1", linewidth=2.5,
            label="Best (final={0:,.0f}, {1} eval.)".format(b_best[-1], len(b_evals)))

    suffix = " ({0})".format(condition_label) if condition_label else ""
    ax.set_title("First vs Best Improvement{0} — seed {1:02d} "
                 "(misma semilla, mismo punto de partida)".format(suffix, seed))
    ax.set_xlabel("iteración (evaluación FMO)")
    ax.set_ylabel("FMO objective")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_dir = os.path.dirname(output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(output, dpi=150)
    plt.close(fig)
    print("seed comparison plot written: {0}".format(output))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--first-csv", required=True, help="trajectory.csv for First Improvement")
    p.add_argument("--best-csv", required=True, help="trajectory.csv for Best Improvement")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--condition-label", default="",
                    help='Appended to the title, e.g. "grid5". Empty for unrestricted.')
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)
    plot_pair(args.first_csv, args.best_csv, args.seed, args.condition_label, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
