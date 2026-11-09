#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the VND (experiments/new/vnd) and general VNS (experiments/new/vns) runs,
with the same figures and tables as the other families (see analyze_new_experiments.py).

  <venv>/bin/python scripts/analyze_vns_experiments.py

Writes figures and CSV tables under experiments/new/{vnd,vns}/comparisons/.
Truncated comparisons cut every run at the same number of real (non-cached) FMO solves.
"""
import os
import sys
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_new_experiments import (NEW, PALETTE, YLAB, load_run, finals, mwu, wilcoxon, fp, f0,
                                     dist, write_csv, grouped_box, conv_panel, save, colormap)

VND = OrderedDict([
    ("VND (desplazamiento 10, 5)", "nangshift10-5"),
    ("VND + swap muestreado", "nangshift10-5-sample"),
    ("VND + swap muestreado, ventana, flip", "nangshift10-5-sample-window-flip"),
])
VNS = OrderedDict([
    ("VNS, perturbación multiescala", "rvns"),
    ("VNS, perturbación de 1 a 3 ángulos", "pershake123"),
    ("VNS con VND (10, 5)", "vnd-rvns"),
    ("VNS con VND + swap muestreado", "vnd-swap"),
])
TABU = OrderedDict([
    ("Tabú, paso 10, first", "tabu/nangshift10/restricted/data/first"),
    ("Tabú, paso 10, best", "tabu/nangshift10/restricted/data/best"),
    ("Tabú, paso 5, first", "tabu/nangshift5/restricted/data/first"),
    ("Tabú, paso 5, best", "tabu/nangshift5/restricted/data/best"),
])
CUTS = (73, 97)
MIN_REL_GAIN, ALPHA = 0.01, 0.05


def runs_of(family, folder):
    base = os.path.join(NEW, family, folder) if folder.startswith("tabu/") else \
        os.path.join(NEW, family, folder, "data", "first")
    if folder.startswith("tabu/"):
        base = os.path.join(NEW, folder)
    return [load_run(os.path.join(base, "seed%02d" % s)) for s in range(1, 16)
            if os.path.exists(os.path.join(base, "seed%02d" % s, "trajectory.csv"))]


def best_at(r, n):
    o = r["objs"][~r["cached"]]
    return float(np.min(o[:n])) if len(o) else float("nan")


def label(gain, p):
    if p < ALPHA and gain >= MIN_REL_GAIN:
        return "positivo"
    if p < ALPHA and gain <= -MIN_REL_GAIN:
        return "negativo"
    return "neutro"


def analyse(family, configs, ref_name, tag):
    out = os.path.join(NEW, family, "comparisons")
    data = OrderedDict((n, runs_of(family, f)) for n, f in configs.items())
    colors = colormap(list(data))
    # convergence (median + IQR), with and without cached rows
    for nc, suf in ((False, ""), (True, "_noncached")):
        fig, ax = plt.subplots(figsize=(8, 5))
        conv_panel(ax, list(data.items()), colors, nc)
        ax.set_title("Convergencia, %s (mediana y rango intercuartil, 15 semillas)" % tag)
        save(fig, os.path.join(out, "V1_conv%s.png" % suf))
    # final objective boxplot
    fig, ax = plt.subplots(figsize=(9, 5.5))
    grouped_box(ax, [(tag, [(n, finals(r), n) for n, r in data.items()])], colors,
                title="Objetivo final por configuración", rotate=20)
    save(fig, os.path.join(out, "V2_final_box.png"))
    # per-run metrics and evaluation counts
    rows = [[n, i + 1, "%.4f" % r["final"], "%.4f" % r["best"], r["n"], r["nn"], r["idx"], r["nn_idx"]]
            for n, rr in data.items() for i, r in enumerate(rr)]
    write_csv(os.path.join(out, "per_run_metrics.csv"),
              ["configuracion", "semilla", "final", "mejor", "evaluaciones", "solves_reales",
               "eval_mejor", "solve_mejor"], rows)
    write_csv(os.path.join(out, "C_eval_counts.csv"),
              ["configuracion", "evals_mediana", "solves_mediana", "solves_min", "solves_max"],
              [[n, f0(np.median([r["n"] for r in rr])), f0(np.median([r["nn"] for r in rr])),
                min(r["nn"] for r in rr), max(r["nn"] for r in rr)] for n, rr in data.items()])
    # statistics against the reference configuration (paired by seed)
    ref = data[ref_name]
    srows = []
    for n, rr in data.items():
        for cut in list(CUTS) + ["final"]:
            a = np.array([r["final"] if cut == "final" else best_at(r, cut) for r in ref])
            b = np.array([r["final"] if cut == "final" else best_at(r, cut) for r in rr])
            k = min(len(a), len(b))
            gain = (a.mean() - b.mean()) / a.mean()
            d = dist(b)
            srows.append([n, cut, k, f0(b.mean()), f0(d["med"]), f0(d["q1"]), f0(d["q3"]),
                          "%.2f" % (100 * gain),
                          "ref" if n == ref_name else fp(wilcoxon(a[:k] - b[:k])),
                          "ref" if n == ref_name else fp(mwu(a, b)),
                          "ref" if n == ref_name else label(gain, wilcoxon(a[:k] - b[:k]))])
    write_csv(os.path.join(out, "V3_stats_vs_%s.csv" % ({"vnd": "vnd_base", "vns": "vns_base", "tabu": "tabu_base"}[family])),
              ["configuracion", "corte_solves", "n", "media", "mediana", "q1", "q3",
               "mejora_pct_vs_ref", "p_wilcoxon_pareado", "p_mann_whitney", "clasificacion"], srows)
    # truncated convergence at a common budget
    cut = CUTS[-1] if family in ("vns", "tabu") else None
    if cut:
        fig, ax = plt.subplots(figsize=(8, 5))
        for n, rr in data.items():
            curves = [np.minimum.accumulate(r["objs"][~r["cached"]])[:cut] for r in rr]
            m = np.array([np.concatenate([c, np.full(cut - len(c), c[-1])]) for c in curves])  # runs that stop early keep their best
            x = np.arange(1, cut + 1)
            ax.plot(x, np.median(m, 0), color=colors[n], label=n, linewidth=1.8)
            ax.fill_between(x, np.percentile(m, 25, 0), np.percentile(m, 75, 0), color=colors[n], alpha=0.2)
        ax.set_xlabel("Evaluaciones no cacheadas (solves FMO reales)")
        ax.set_ylabel("Mejor objetivo hasta el momento (FMO)")
        ax.set_title("Convergencia truncada a %d solves (mediana y rango intercuartil)" % cut)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        save(fig, os.path.join(out, "V1_conv_truncada.png"))
    print(family, {n: len(r) for n, r in data.items()})


if __name__ == "__main__":
    analyse("vnd", VND, "VND (desplazamiento 10, 5)", "VND")
    analyse("vns", VNS, "VNS, perturbación multiescala", "VNS")
    analyse("tabu", TABU, "Tabú, paso 10, first", "Tabú")
