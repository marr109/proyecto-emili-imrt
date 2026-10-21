#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grouped boxplot: N boxes organized into GROUPS (e.g. catalog: restricted /
unrestricted), each box colored consistently by VARIANT (e.g. fijo / random /
circular) across the whole figure, with a legend mapping color -> variant
instead of forcing the reader to decode long text labels per box.

Each --test is (LABEL, BASE_DIR, STRATEGY, GROUP, VARIANT):
  LABEL    short x-tick text (e.g. "First", "Best")
  BASE_DIR directory containing data/<STRATEGY>/seedNN/trajectory.csv
  STRATEGY "first" or "best"
  GROUP    cluster this box belongs to (e.g. "Restricted", "Unrestricted") --
           boxes are ordered by group first, a visual gap separates groups,
           and each group gets a centered label below the x-axis.
  VARIANT  drives color + legend entry (e.g. "Fijo", "Random", "Circular") --
           same VARIANT string always gets the same color across the figure.

CLI example:
  python scripts/plot_boxplot_grouped.py \
      --test "First" DIR1 first Restricted Fijo \
      --test "Best"  DIR1 best  Restricted Fijo \
      --seeds 1-15 --title "..." --output out.png
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
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory

# Fixed categorical order (validated colorblind-safe, matplotlib tab10) --
# assigned by VARIANT in first-seen order, never cycled per-box.
CATEGORICAL_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

GROUP_GAP = 1.0  # extra x-space inserted between consecutive groups


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
    p.add_argument("--test", action="append", nargs=5, required=True,
                    metavar=("LABEL", "BASE_DIR", "STRATEGY", "GROUP", "VARIANT"),
                    help="Repeatable: one box per --test.")
    p.add_argument("--seeds", required=True, help='e.g. "1-15" or "1,2,3,11-15"')
    p.add_argument("--title", required=True)
    p.add_argument("--ylabel", default="Objetivo final (FMO)")
    p.add_argument("--figsize", default=None,
                   help="ancho,alto en pulgadas; por defecto se deduce del numero de cajas")
    p.add_argument("--legend-title", default="Orden de escaneo",
                   help="encabezado del recuadro de leyenda")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    seeds = parse_seeds(args.seeds)

    # Assign color per VARIANT in first-seen order (stable across the figure).
    variant_order = []
    for _, _, _, _, variant in args.test:
        if variant not in variant_order:
            variant_order.append(variant)
    if len(variant_order) > len(CATEGORICAL_COLORS):
        raise SystemExit(
            "{0} variants but only {1} colors in the fixed palette -- extend "
            "CATEGORICAL_COLORS deliberately, don't cycle it".format(
                len(variant_order), len(CATEGORICAL_COLORS)))
    variant_color = {v: CATEGORICAL_COLORS[i] for i, v in enumerate(variant_order)}

    # Preserve GROUP order as first-seen, keep --test order within a group.
    group_order = []
    for _, _, _, group, _ in args.test:
        if group not in group_order:
            group_order.append(group)

    tests_by_group = {g: [] for g in group_order}
    for label, base, strategy, group, variant in args.test:
        tests_by_group[group].append((label, base, strategy, variant))

    labels, datasets, colors = [], [], []
    xpos = []
    group_centers = {}
    x = 1.0
    for group in group_order:
        start_x = x
        for label, base, strategy, variant in tests_by_group[group]:
            values = [final_objective(base, strategy, s) for s in seeds]
            labels.append(label)
            datasets.append(values)
            colors.append(variant_color[variant])
            xpos.append(x)
            x += 1.0
        group_centers[group] = (start_x + (x - 1.0)) / 2.0
        x += GROUP_GAP

    n_boxes = len(labels)
    if args.figsize:
        fw, fh = [float(v) for v in args.figsize.split(",")]
    else:
        fw, fh = max(7.0, 1.1 * n_boxes + 2.5), 6.5
    fig, ax = plt.subplots(figsize=(fw, fh))

    bp = ax.boxplot(
        datasets, positions=xpos, **{_BOXPLOT_LABEL_KW: labels}, patch_artist=True, widths=0.6,
        medianprops=dict(color="#2b2b2b", linewidth=1.6),
        whiskerprops=dict(color="#555555", linewidth=1.1),
        capprops=dict(color="#555555", linewidth=1.1),
        boxprops=dict(linewidth=1.1),
        flierprops=dict(marker="o", markersize=4, markerfacecolor="none",
                         markeredgecolor="#555555", alpha=0.7),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)

    # Overlay individual seed points (jittered) -- n is small (~15), showing
    # the raw points alongside the summary avoids hiding structure the box
    # alone would smooth over.
    rnd = random.Random(12345)
    for xi, values, color in zip(xpos, datasets, colors):
        xs = [xi + rnd.uniform(-0.14, 0.14) for _ in values]
        ax.scatter(xs, values, color=color, s=22, alpha=0.8, zorder=3,
                   edgecolors="white", linewidths=0.4)

    # Vertical separators + centered group labels below the x-tick labels.
    ax.set_xlim(0.3, x - GROUP_GAP + 0.3)
    y0, y1 = ax.get_ylim()
    if len(group_order) > 1:
        boundary_x = None
        prev_end = None
        for group in group_order:
            center = group_centers[group]
            if prev_end is not None:
                boundary_x = (prev_end + (center - (len(tests_by_group[group]) - 1) / 2.0)) / 2.0
                ax.axvline(boundary_x, color="#cccccc", linewidth=1.0, zorder=0)
            prev_end = center + (len(tests_by_group[group]) - 1) / 2.0

    for group, center in group_centers.items():
        ax.text(center, -0.16, group, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=11, fontweight="bold")

    legend_handles = [Patch(facecolor=variant_color[v], edgecolor=variant_color[v],
                             alpha=0.55, label=v) for v in variant_order]
    ax.legend(handles=legend_handles, title=args.legend_title, loc="upper right",
              frameon=True, framealpha=0.9, fontsize=9, title_fontsize=9)

    ax.set_ylabel(args.ylabel)
    ax.set_title(args.title, wrap=True)
    ax.grid(alpha=0.25, axis="y")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    ax.text(0.99, 0.02, "n={0} semillas por caja".format(len(seeds)),
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="#777777")

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.16)
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)

    print("boxplot written: {0}".format(args.output))
    for label, group, variant, values in zip(
            labels, [g for g in group_order for _ in tests_by_group[g]], [v for g in group_order for (_, _, _, v) in tests_by_group[g]], datasets):
        s = sorted(values)
        median = s[len(s) // 2] if len(s) % 2 else (s[len(s)//2 - 1] + s[len(s)//2]) / 2
        print("  {0}/{1}/{2}: n={3} median={4:.1f} mean={5:.1f} min={6:.1f} max={7:.1f}".format(
            group, label, variant, len(values), median, sum(values) / len(values), min(values), max(values)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
