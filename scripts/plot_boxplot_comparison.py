#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generic boxplot comparing N "tests" (local-search configurations run across
multiple seeds), showing the distribution of the final FMO objective per
test. Built to take an arbitrary number of tests — two today (First vs
Best), ten tomorrow — not hardcoded to any fixed comparison.

Each --test is (LABEL, BASE_DIR, STRATEGY); the script reads
<BASE_DIR>/data/<STRATEGY>/seedNN/trajectory.csv for every seed in --seeds
and takes the minimum objective (the local optimum reached) as that seed's
data point.

CLI example:
  python scripts/plot_boxplot_comparison.py \
      --test "First (grid5)" experiments/local_search/nangshift10/restricted first \
      --test "Best (grid5)"  experiments/local_search/nangshift10/restricted best \
      --seeds 1-15 \
      --title "grid5, nangshift 10 — distribución del objetivo final" \
      --output experiments/local_search/nangshift10/restricted/plots/boxplot.png
"""

import argparse
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# matplotlib renombro 'labels' -> 'tick_labels' en 3.9. calafate corre 3.3
# (ligada a Python 3.6), asi que el nombre se elige en runtime.
_BOXPLOT_LABEL_KW = ("tick_labels"
                     if tuple(int(x) for x in matplotlib.__version__.split(".")[:2]) >= (3, 9)
                     else "labels")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory

# Fixed categorical order (validated colorblind-safe, matplotlib tab10) —
# assigned by input order, never cycled. Extend deliberately if a
# comparison ever needs more than 10 boxes.
CATEGORICAL_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def parse_seeds(spec):
    seeds = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            seeds.extend(range(int(a), int(b) + 1))
        else:
            seeds.append(int(part))
    return seeds


def final_objective(base, strategy, seed):
    path = os.path.join(base, "data", strategy, "seed{0:02d}".format(seed), "trajectory.csv")
    _, objectives, _ = read_trajectory(path)
    return min(objectives)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--test", action="append", nargs=3, required=True,
                    metavar=("LABEL", "BASE_DIR", "STRATEGY"),
                    help="Repeatable: one box per --test.")
    p.add_argument("--seeds", required=True, help='e.g. "1-15" or "1,2,3,11-15"')
    p.add_argument("--title", required=True)
    p.add_argument("--ylabel", default="FMO objective final")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    seeds = parse_seeds(args.seeds)
    if len(args.test) > len(CATEGORICAL_COLORS):
        raise SystemExit(
            "{0} tests but only {1} colors in the fixed palette -- extend "
            "CATEGORICAL_COLORS deliberately, don't cycle it".format(
                len(args.test), len(CATEGORICAL_COLORS)))

    labels, datasets = [], []
    for label, base, strategy in args.test:
        values = [final_objective(base, strategy, s) for s in seeds]
        labels.append(label)
        datasets.append(values)

    fig_width = max(6.0, 1.4 * len(labels) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 6.5))

    bp = ax.boxplot(
        datasets, **{_BOXPLOT_LABEL_KW: labels}, patch_artist=True, widths=0.5,
        medianprops=dict(color="#2b2b2b", linewidth=1.6),
        whiskerprops=dict(color="#555555", linewidth=1.1),
        capprops=dict(color="#555555", linewidth=1.1),
        boxprops=dict(linewidth=1.1),
        flierprops=dict(marker="o", markersize=4, markerfacecolor="none",
                         markeredgecolor="#555555", alpha=0.7),
    )
    for patch, color in zip(bp["boxes"], CATEGORICAL_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)

    # Overlay individual seed points (jittered) -- n is small (≈15), showing
    # the raw points alongside the summary avoids hiding structure the box
    # alone would smooth over.
    rnd = random.Random(12345)
    for i, (values, color) in enumerate(zip(datasets, CATEGORICAL_COLORS), start=1):
        xs = [i + rnd.uniform(-0.12, 0.12) for _ in values]
        ax.scatter(xs, values, color=color, s=22, alpha=0.8, zorder=3,
                   edgecolors="white", linewidths=0.4)

    ax.set_ylabel(args.ylabel)
    ax.set_title(args.title, wrap=True)
    ax.grid(alpha=0.25, axis="y")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    ax.text(0.99, 0.02, "n={0} semillas por caja".format(len(seeds)),
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="#777777")

    fig.tight_layout()
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)

    print("boxplot written: {0}".format(args.output))
    for label, values in zip(labels, datasets):
        s = sorted(values)
        median = s[len(s) // 2] if len(s) % 2 else (s[len(s)//2 - 1] + s[len(s)//2]) / 2
        print("  {0}: n={1} median={2:.1f} mean={3:.1f} min={4:.1f} max={5:.1f}".format(
            label, len(values), median, sum(values) / len(values), min(values), max(values)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
