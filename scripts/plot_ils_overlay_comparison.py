#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overlay de las N trayectorias (mejor-hasta-el-momento) de First vs Best,
un panel por estrategia, para ver dispersión entre semillas en vez de solo
el punto final (que ya cubre el boxplot/dumbbell).

En el panel de First se marcan además las perturbaciones detectadas
(mismo criterio que plot_ils_trajectory.py). Best no se marca: la búsqueda
local Best no rebasa a mitad de ronda, así que el detector (pensado para
First) no aplica — ver docstring de detect_perturbations.

CLI example:
  python scripts/plot_ils_overlay_comparison.py \
      --base experiments/ils/nangshift5/grid5 \
      --seeds 1-15 --step 5 \
      --title "ILS First vs Best — nangshift 5, grid5" \
      --output experiments/ils/nangshift5/grid5/plots/overlay_first_vs_best.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory, running_best
from plot_ils_trajectory import detect_perturbations

# Mismo orden categórico fijo que plot_boxplot_comparison.py (tab10,
# validado colorblind-safe) — First=azul, Best=naranja, nunca ciclado.
COLOR_FIRST = "#1f77b4"
COLOR_BEST = "#ff7f0e"
COLOR_PERTURBATION = "#d1495b"


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


def plot_panel(ax, base, strategy, seeds, color, step, mark_perturbations, title):
    n_plotted = 0
    n_perturbations = 0
    for seed in seeds:
        path = os.path.join(base, "data", strategy, "seed{0:02d}".format(seed), "trajectory.csv")
        if not os.path.exists(path):
            continue
        evals, objectives, angle_sets = read_trajectory(path)
        best = running_best(objectives)
        ax.step(evals, best, where="post", color=color, linewidth=1.1, alpha=0.35, zorder=2)
        n_plotted += 1

        if mark_perturbations:
            records = detect_perturbations(evals, objectives, angle_sets, step, reject_repeated=True)
            pert_set = {r["eval_despues"] for r in records}
            pert_evals = [e for e in evals if e in pert_set]
            pert_best = [b for e, b in zip(evals, best) if e in pert_set]
            ax.scatter(pert_evals, pert_best, color=COLOR_PERTURBATION, s=10,
                       alpha=0.6, zorder=3, linewidths=0)
            n_perturbations += len(pert_evals)

    ax.set_title("{0} (n={1} semillas)".format(title, n_plotted))
    ax.set_xlabel("evaluación")
    ax.grid(alpha=0.3)
    return n_plotted, n_perturbations


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base", required=True, help="<base>/data/{first,best}/seedNN/trajectory.csv")
    p.add_argument("--seeds", required=True, help='e.g. "1-15" or "1,2,3,11-15"')
    p.add_argument("--step", type=int, required=True, help="step de nangshift usado (5 o 10)")
    p.add_argument("--title", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    seeds = parse_seeds(args.seeds)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    n1, npert1 = plot_panel(ax1, args.base, "first", seeds, COLOR_FIRST, args.step,
                             mark_perturbations=True, title="First")
    n2, _ = plot_panel(ax2, args.base, "best", seeds, COLOR_BEST, args.step,
                        mark_perturbations=False, title="Best")

    ax1.set_ylabel("FMO objective (mejor hasta el momento)")

    line_first = plt.Line2D([0], [0], color=COLOR_FIRST, linewidth=1.5, alpha=0.7, label="trayectoria por semilla")
    dot_pert = plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COLOR_PERTURBATION,
                          markersize=6, label="perturbación ({0})".format(npert1))
    ax1.legend(handles=[line_first, dot_pert], loc="upper right")

    line_best = plt.Line2D([0], [0], color=COLOR_BEST, linewidth=1.5, alpha=0.7, label="trayectoria por semilla")
    ax2.legend(handles=[line_best], loc="upper right")
    ax2.set_title(ax2.get_title() + "  (perturbaciones no marcadas — detector solo válido para First)",
                  fontsize=9)

    fig.suptitle(args.title)
    fig.tight_layout()
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print("overlay comparison written: {0} (First n={1}, Best n={2}, perturbaciones First={3})".format(
        args.output, n1, n2, npert1))


if __name__ == "__main__":
    sys.exit(main())
