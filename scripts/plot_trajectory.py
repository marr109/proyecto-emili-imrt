#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the BAO search trajectory from a CSV produced by `BaoProblem::openCsvLog`
(columns: eval,angles_deg,objective).

Two panels:
  1. Objective value per evaluation, with a running-best (incumbent) overlay.
  2. Heatmap of which catalog angles were active at each evaluation.

CLI example:
  python scripts/plot_trajectory.py --csv /tmp/bao_shift_trajectory.csv \
      --output plots/trajectory.png \
      --improvements-csv trajectory_improvements.csv \
      --convergence-output plots/convergence.png
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_trajectory(path):
    evals, objectives, angle_sets = [], [], []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            evals.append(int(row["eval"]))
            objectives.append(float(row["objective"]))
            angle_sets.append([int(a) for a in row["angles_deg"].split(";")])
    return evals, objectives, angle_sets


def running_best(objectives):
    best = []
    cur = float("inf")
    for v in objectives:
        cur = min(cur, v)
        best.append(cur)
    return best


def improving_records(evals, objectives, angle_sets):
    """Return the initial solution and every subsequent strict record improvement."""
    records = []
    best = float("inf")
    for ev, objective, angles in zip(evals, objectives, angle_sets):
        if objective < best:
            records.append((ev, list(angles), objective))
            best = objective
    return records


def write_improvements(path, records):
    out_dir = os.path.dirname(path) or "."
    os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["improvement", "eval", "angles_deg", "objective"])
        for number, (ev, angles, objective) in enumerate(records):
            writer.writerow([
                number,
                ev,
                ";".join(str(a) for a in angles),
                "{0:.6f}".format(objective),
            ])


def plot_convergence(path, records):
    out_dir = os.path.dirname(path) or "."
    os.makedirs(out_dir, exist_ok=True)

    improvement_numbers = list(range(len(records)))
    objectives = [record[2] for record in records]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(
        improvement_numbers, objectives,
        color="#2f6f9f", marker="o", linewidth=2,
    )
    for number, objective in zip(improvement_numbers, objectives):
        ax.annotate(
            "{0:.0f}".format(objective),
            (number, objective), xytext=(0, 7),
            textcoords="offset points", ha="center", fontsize=8,
        )
    ax.set_xlabel("mejora aceptada")
    ax.set_ylabel("FMO objective")
    ax.set_title("Convergencia BAO — soluciones que mejoran el incumbente")
    ax.set_xticks(improvement_numbers)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_iterations(path, evals, objectives):
    out_dir = os.path.dirname(path) or "."
    os.makedirs(out_dir, exist_ok=True)

    best = running_best(objectives)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(evals, objectives, color="#9aa5b1", linewidth=0.8, alpha=0.6, label="evaluado")
    ax.step(evals, best, where="post", color="#2f6f9f", linewidth=2, label="mejor hasta el momento")
    ax.set_xlabel("iteración (evaluación FMO)")
    ax.set_ylabel("FMO objective")
    ax.set_title("Objetivo por iteración")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def build_activity_matrix(angle_sets):
    catalog = sorted({a for s in angle_sets for a in s})
    index = {a: i for i, a in enumerate(catalog)}
    mat = np.zeros((len(catalog), len(angle_sets)))
    for j, s in enumerate(angle_sets):
        for a in s:
            mat[index[a], j] = 1.0
    return catalog, mat


def main(argv=None):
    p = argparse.ArgumentParser(description="Plot a BAO search trajectory.")
    p.add_argument("--csv", required=True, help="Trajectory CSV (eval,angles_deg,objective).")
    p.add_argument("--output", default="plots/trajectory.png", help="Output PNG path.")
    p.add_argument(
        "--skip-trajectory",
        action="store_true",
        help="Skip the two-panel trajectory.png (objective+heatmap by evaluation); "
             "only improvements-csv/convergence-output are produced.",
    )
    p.add_argument(
        "--improvements-csv",
        help="Optional CSV containing only the initial solution and strict improvements.",
    )
    p.add_argument(
        "--convergence-output",
        help="Optional convergence PNG using only strict improvements.",
    )
    p.add_argument(
        "--iterations-output",
        help="Optional PNG: objective (+ running best) per iteration/evaluation, no heatmap.",
    )
    args = p.parse_args(argv)

    evals, objectives, angle_sets = read_trajectory(args.csv)
    if not evals:
        raise SystemExit("error: empty trajectory CSV: {0}".format(args.csv))

    best = running_best(objectives)
    improvements = improving_records(evals, objectives, angle_sets)
    catalog, activity = build_activity_matrix(angle_sets)

    if args.improvements_csv:
        write_improvements(args.improvements_csv, improvements)
    if args.convergence_output:
        plot_convergence(args.convergence_output, improvements)
    if args.iterations_output:
        plot_iterations(args.iterations_output, evals, objectives)

    if args.skip_trajectory:
        if args.improvements_csv:
            print("improvements CSV written: {0}".format(args.improvements_csv))
        if args.convergence_output:
            print("convergence plot written: {0}".format(args.convergence_output))
        if args.iterations_output:
            print("iterations plot written: {0}".format(args.iterations_output))
        print("  evaluations : {0}".format(len(evals)))
        print("  improvements: {0}".format(len(improvements) - 1))
        print("  best objective: {0:.2f} (eval {1})".format(
            min(objectives), evals[objectives.index(min(objectives))]))
        return 0

    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1.4]},
    )

    ax1.plot(evals, objectives, color="#9aa5b1", linewidth=0.8, alpha=0.6, label="evaluado")
    ax1.step(evals, best, where="post", color="#2f6f9f", linewidth=2, label="mejor hasta el momento")
    ax1.set_ylabel("FMO objective")
    ax1.set_title("Trayectoria de búsqueda BAO")
    ax1.legend(loc="upper right")
    ax1.grid(alpha=0.3)

    ax2.imshow(
        activity, aspect="auto", cmap="Greys", interpolation="nearest",
        extent=[evals[0] - 0.5, evals[-1] + 0.5, -0.5, len(catalog) - 0.5],
        origin="lower",
    )
    ax2.set_yticks(range(len(catalog)))
    ax2.set_yticklabels([str(a) for a in catalog], fontsize=7)
    ax2.set_ylabel("ángulo activo (°)")
    ax2.set_xlabel("evaluación")

    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print("trajectory plot written: {0}".format(args.output))
    if args.improvements_csv:
        print("improvements CSV written: {0}".format(args.improvements_csv))
    if args.convergence_output:
        print("convergence plot written: {0}".format(args.convergence_output))
    if args.iterations_output:
        print("iterations plot written: {0}".format(args.iterations_output))
    print("  evaluations : {0}".format(len(evals)))
    print("  improvements: {0}".format(len(improvements) - 1))
    print("  best objective: {0:.2f} (eval {1})".format(
        min(objectives), evals[objectives.index(min(objectives))]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
