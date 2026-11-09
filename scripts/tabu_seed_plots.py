#!/usr/bin/env python3
"""Write convergence.png and trajectory.png next to each trajectory.csv of the tabu and tabu_attr runs.

convergence.png: best objective so far vs real FMO solves (from convergence.csv), with every evaluated candidate as dots.
trajectory.png:  objective of the current solution at each tabu iteration (accepted moves, from tabu_events.csv),
                 the candidates rejected as tabu, and the best objective so far.

Needs convergence.csv and tabu_events.csv (tabu_seed_convergence_csv.py, tabu_replay_events.py).

Usage: tabu_seed_plots.py [experiments/new]
"""
import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def convergence(d, root):
    rows = list(csv.DictReader(open(os.path.join(d, "convergence.csv"))))
    if not rows:
        return
    xs = [int(r["real_solves"]) for r in rows if r["cached"] != "true" and float(r["objective"]) < 1e20]
    ys = [float(r["objective"]) for r in rows if r["cached"] != "true" and float(r["objective"]) < 1e20]
    bx = [max(int(r["real_solves"]), 1) for r in rows]
    by = [float(r["best_so_far"]) for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(xs, ys, ".", color="#999999", ms=3, label="objetivo de cada candidato")
    ax.step(bx, by, where="post", color="black", lw=1.6, label="mejor objetivo hasta el momento (%.0f)" % by[-1])
    ax.set_xlabel("Resoluciones del FMO")
    ax.set_ylabel("Objetivo del FMO")
    ax.grid(alpha=.3, lw=.4)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(os.path.relpath(d, root), fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "convergence.png"), dpi=150)
    plt.close(fig)


def trajectory(d, root):
    ev = list(csv.DictReader(open(os.path.join(d, "tabu_events.csv"))))
    mv = [(int(e["iter"]), float(e["objective"])) for e in ev if e["event"] in ("start", "move")]
    sk = [(int(e["iter"]), float(e["objective"])) for e in ev if e["event"] == "tabu_skip" and float(e["objective"]) < 1e20]
    best, bx, by = float("inf"), [], []
    for it, o in mv:
        best = min(best, o)
        bx.append(it)
        by.append(best)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    if sk:
        ax.plot([i for i, _ in sk], [o for _, o in sk], "x", color="#CC79A7", ms=4, label="candidato rechazado por tabú")
    ax.plot([i for i, _ in mv], [o for _, o in mv], "-o", color="#0072B2", ms=3, lw=1.2, label="solución actual")
    ax.step(bx, by, where="post", color="black", lw=1.4, label="mejor objetivo hasta el momento (%.0f)" % best)
    ax.set_xlabel("Iteración de la búsqueda tabú")
    ax.set_ylabel("Objetivo del FMO")
    ax.grid(alpha=.3, lw=.4)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(os.path.relpath(d, root), fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "trajectory.png"), dpi=150)
    plt.close(fig)


def main(root):
    n = 0
    for d, _, files in sorted(os.walk(root)):
        if "convergence.csv" not in files or "tabu_events.csv" not in files:
            continue
        if not re.search(r"/(tabu_attr|tabu)/nangshift", d + "/") or "catalog_restricted" in d or "tabu_attr_trace" in d:
            continue
        convergence(d, root)
        trajectory(d, root)
        n += 1
    print("seeds plotted:", n)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "experiments/new")
