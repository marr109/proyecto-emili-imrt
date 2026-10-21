#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aggregate First vs Best Improvement comparison across seeds:
  1. Final objective per seed (dumbbell plot).
  2. Evaluation cost per seed (grouped bar chart).

Two modes:
  --mode strategy-pair
      One condition, First vs Best per seed (2 series).
      Reads <base>/data/{first,best}/seedNN/trajectory.csv

  --mode condition-pair
      Two conditions (e.g. unrestricted vs grid5), each with First and Best
      (4 series). Connects same-strategy points across conditions to show
      the effect of the condition change.
      Reads <condition-a-base>/data/{first,best}/seedNN/trajectory.csv and
            <condition-b-base>/data/{first,best}/seedNN/trajectory.csv

CLI examples:
  python scripts/plot_aggregate_comparison.py --mode strategy-pair \
      --base experiments/local_search/nangshift10/unrestricted \
      --seeds 1-15 \
      --title-top "First vs Best Improvement — objetivo final por semilla (misma semilla en ambos, locmin, nangshift 10)" \
      --title-bottom "Costo de cómputo — cantidad de evaluaciones hasta óptimo local" \
      --output experiments/local_search/nangshift10/unrestricted/plots/first_vs_best_comparison.png

  python scripts/plot_aggregate_comparison.py --mode condition-pair \
      --condition-a-label "unrestricted (360)" --condition-a-base experiments/local_search/nangshift10/unrestricted \
      --condition-b-label "grid5 (múltiplos de 5°)" --condition-b-base experiments/local_search/nangshift10/restricted \
      --seeds 1-15 \
      --title-top "Objetivo final por semilla — unrestricted (360°) vs grid5 (múltiplos de 5°)" \
      --title-bottom "Costo de cómputo hasta óptimo local" \
      --output experiments/local_search/nangshift10/comparisons/unrestricted_vs_grid5/plots/grid5_vs_unrestricted_comparison.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory


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


def final_and_evals(csv_path):
    evals, objectives, _ = read_trajectory(csv_path)
    return min(objectives), len(evals)


def seed_csv(base, strategy, seed):
    return os.path.join(base, "data", strategy, "seed{0:02d}".format(seed), "trajectory.csv")


def strategy_pair_mode(args):
    seeds = args.seeds
    first_obj, first_ev, best_obj, best_ev = [], [], [], []
    for s in seeds:
        o, e = final_and_evals(seed_csv(args.base, "first", s))
        first_obj.append(o)
        first_ev.append(e)
        o, e = final_and_evals(seed_csv(args.base, "best", s))
        best_obj.append(o)
        best_ev.append(e)

    labels = ["seed {0:02d}".format(s) for s in seeds]
    x = list(range(len(seeds)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(10, 1.1 * len(seeds)), 9))

    for xi, fo, bo in zip(x, first_obj, best_obj):
        ax1.plot([xi, xi], [fo, bo], color="gray", linewidth=1, zorder=1)
    ax1.scatter(x, first_obj, color="C0", s=70, label="First Improvement", zorder=2)
    ax1.scatter(x, best_obj, color="C1", s=70, label="Best Improvement", zorder=2)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("Mejor objetivo FMO alcanzado")
    ax1.set_title(args.title_top)
    ax1.legend(loc="best")
    ax1.grid(alpha=0.3, axis="y")

    width = 0.35
    ax2.bar([xi - width / 2 for xi in x], first_ev, width, color="C0", label="First Improvement")
    ax2.bar([xi + width / 2 for xi in x], best_ev, width, color="C1", label="Best Improvement")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Evaluaciones FMO (costo)")
    ax2.set_title(args.title_bottom)
    ax2.grid(alpha=0.3, axis="y")

    _save(fig, args.output)


def condition_pair_mode(args):
    seeds = args.seeds
    a_label, b_label = args.condition_a_label, args.condition_b_label
    data = {}
    for cond_label, base in ((a_label, args.condition_a_base), (b_label, args.condition_b_base)):
        for strat in ("first", "best"):
            objs, evs = [], []
            for s in seeds:
                o, e = final_and_evals(seed_csv(base, strat, s))
                objs.append(o)
                evs.append(e)
            data[(cond_label, strat)] = (objs, evs)

    labels = ["seed {0:02d}".format(s) for s in seeds]
    n = len(seeds)
    x = list(range(n))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(10, 1.3 * n), 10))

    jitter = 0.12
    for strat, color in (("first", "C0"), ("best", "C1")):
        a_objs, _ = data[(a_label, strat)]
        b_objs, _ = data[(b_label, strat)]
        xoff = -jitter if strat == "first" else jitter
        xs = [xi + xoff for xi in x]
        for xi, ao, bo in zip(xs, a_objs, b_objs):
            ax1.plot([xi, xi], [ao, bo], color=color, alpha=0.5, linewidth=1, zorder=1)
        ax1.scatter(xs, a_objs, facecolors="none", edgecolors=color, s=70, zorder=2,
                    label="{0} — {1}".format(strat.capitalize(), a_label))
        ax1.scatter(xs, b_objs, facecolors=color, edgecolors=color, s=70, zorder=2,
                    label="{0} — {1}".format(strat.capitalize(), b_label))
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("FMO objective final")
    ax1.set_title(args.title_top)
    ax1.legend(loc="best", ncol=2, fontsize=9)
    ax1.grid(alpha=0.3, axis="y")

    width = 0.19
    slots = {("first", a_label): -1.5, ("first", b_label): -0.5,
             ("best", a_label): 0.5, ("best", b_label): 1.5}
    for strat, color in (("first", "C0"), ("best", "C1")):
        for cond_label in (a_label, b_label):
            _, evs = data[(cond_label, strat)]
            off = slots[(strat, cond_label)]
            xs = [xi + off * width for xi in x]
            if cond_label == a_label:
                ax2.bar(xs, evs, width, fill=False, edgecolor=color, hatch="//",
                        label="{0} — {1}".format(strat.capitalize(), cond_label))
            else:
                ax2.bar(xs, evs, width, color=color,
                        label="{0} — {1}".format(strat.capitalize(), cond_label))
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Evaluaciones FMO (costo)")
    ax2.set_title(args.title_bottom)
    ax2.legend(loc="best", ncol=2, fontsize=9)
    ax2.grid(alpha=0.3, axis="y")

    _save(fig, args.output)


def _save(fig, output):
    fig.tight_layout()
    out_dir = os.path.dirname(output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(output, dpi=150)
    plt.close(fig)
    print("aggregate plot written: {0}".format(output))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["strategy-pair", "condition-pair"], required=True)
    p.add_argument("--seeds", required=True, help='e.g. "1-15" or "1,2,3,11-15"')
    p.add_argument("--title-top", required=True)
    p.add_argument("--title-bottom", required=True)
    p.add_argument("--output", required=True)
    # strategy-pair
    p.add_argument("--base", help="experiment base dir (strategy-pair mode)")
    # condition-pair
    p.add_argument("--condition-a-label")
    p.add_argument("--condition-a-base")
    p.add_argument("--condition-b-label")
    p.add_argument("--condition-b-base")
    args = p.parse_args(argv)
    args.seeds = parse_seeds(args.seeds)

    if args.mode == "strategy-pair":
        if not args.base:
            raise SystemExit("--base is required for --mode strategy-pair")
        strategy_pair_mode(args)
    else:
        missing = [n for n in ("condition_a_label", "condition_a_base",
                                 "condition_b_label", "condition_b_base")
                   if not getattr(args, n)]
        if missing:
            raise SystemExit("--mode condition-pair requires: " + ", ".join(
                "--" + m.replace("_", "-") for m in missing))
        condition_pair_mode(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
