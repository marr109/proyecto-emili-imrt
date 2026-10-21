#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of experiments/new (local search, SA, ILS; per-organ FMO model, K=4).

Run with a Python that has numpy, scipy and matplotlib, e.g.
  <venv>/bin/python scripts/analyze_new_experiments.py

Writes figures and CSV tables under experiments/new/<family>/comparisons/ and
the summary experiments/new/ANALYSIS_NOTES.md. Only reads trajectory.csv and
run.log of the runs; it does not read any plans_dosimetry.csv.
"""
import csv
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEW = os.path.join(ROOT, "experiments", "new")
STEPS = [5, 10]
CATS = ["restricted", "unrestricted"]
CAT_ES = {"restricted": "Restricted", "unrestricted": "Unrestricted"}
SEEDS = list(range(1, 16))
ORDERS = OrderedDict([("fixed", ""), ("circular", "-shift-circular"), ("random", "-shift-rand")])
ORDER_ES = {"fixed": "Fijo", "circular": "Circular", "random": "Random"}
PERTS = OrderedDict([("range3", ""), ("adaptive", "-adaptive"), ("mag2", "-mag2"),
                     ("mag2-stepwise", "-mag2-stepwise"),
                     ("adaptivemag-stepwise", "-adaptivemag-stepwise")])
PERT_ES = {"range3": "range3 (3 ang.)", "adaptive": "adaptive", "mag2": "mag2 (2 ang.)",
           "mag2-stepwise": "mag2-stepwise", "adaptivemag-stepwise": "adaptivemag-stepwise"}
CALS = OrderedDict([("linear", "-budget"), ("geometric", "-budget-geometric")])
CAL_ES = {"linear": "SA lineal", "geometric": "SA geométrico"}
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
YLAB = "Objetivo final (FMO)"

_cache = {}


# ----------------------------------------------------------------- loading
def load_run(d):
    if d in _cache:
        return _cache[d]
    objs, cached = [], []
    with open(os.path.join(d, "trajectory.csv")) as f:
        for r in csv.DictReader(f):
            objs.append(float(r["objective"]))
            cached.append(r["cached"].strip().lower() == "true")
    objs = np.array(objs)
    cached = np.array(cached)
    final = None
    try:
        with open(os.path.join(d, "run.log")) as f:
            for line in f:
                m = re.search(r"FMO objective f\*\s*:\s*([-+0-9.eE]+)", line)
                if m:
                    final = float(m.group(1))
    except IOError:
        pass
    best = float(objs.min())
    idx = int(np.argmax(objs == best))  # 0-based first index of the min
    run = dict(objs=objs, cached=cached, n=len(objs), nc=int(cached.sum()),
               nn=int((~cached).sum()), best=best, final=final if final is not None else best,
               has_final=final is not None, idx=idx + 1,
               nn_idx=int((~cached[:idx + 1]).sum()), first_obj=float(objs[0]))
    _cache[d] = run
    return run


def ls_dir(step, cat, order, strat, seed):
    return os.path.join(NEW, "local_search", "nangshift%d" % step, cat + ORDERS[order], "data",
                        strat, "seed%02d" % seed)


def sa_dir(step, cat, cal, seed):
    return os.path.join(NEW, "sa", "nangshift%d" % step, cat + CALS[cal], "data", "seed%02d" % seed)


def ils_dir(step, cat, pert, order, strat, seed):
    return os.path.join(NEW, "ils", "nangshift%d" % step, cat + PERTS[pert] + ORDERS[order], "data",
                        strat, "seed%02d" % seed)


def ls_runs(step, cat, order, strat):
    return [load_run(ls_dir(step, cat, order, strat, s)) for s in SEEDS]


def sa_runs(step, cat, cal):
    return [load_run(sa_dir(step, cat, cal, s)) for s in SEEDS]


def ils_runs(step, cat, pert, order, strat):
    return [load_run(ils_dir(step, cat, pert, order, strat, s)) for s in SEEDS]


def finals(runs):
    return np.array([r["final"] for r in runs])


def anytime(runs):
    return np.array([r["best"] for r in runs])


# ------------------------------------------------------------------- stats
def fp(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n/d"
    return "%.4f" % p if p >= 0.0001 else "<0.0001"


def mwu(a, b):
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def kw(groups):
    try:
        return float(stats.kruskal(*groups).pvalue)
    except ValueError:
        return float("nan")


def wilcoxon(d):
    d = np.asarray(d)
    if np.all(d == 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(d).pvalue)
    except ValueError:
        return float("nan")


def q(a, p):
    return float(np.percentile(a, p))


def dist(a):
    a = np.asarray(a, dtype=float)
    return dict(med=q(a, 50), q1=q(a, 25), q3=q(a, 75), p90=q(a, 90), mn=a.min(), mx=a.max())


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out) + "\n"


def write_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def f0(x):
    return "%.0f" % x


def f1(x):
    return "%.1f" % x


# ---------------------------------------------------------------- plotting
def grouped_box(ax, groups, colors, ylabel=YLAB, title=None, legend_title=None, rotate=0):
    """groups: list of (group_label, [(box_label, values, variant)])."""
    x = 1.0
    xs, labels, data, cols, centers = [], [], [], [], []
    variants = []
    for gname, boxes in groups:
        start = x
        for lab, vals, var in boxes:
            xs.append(x)
            labels.append(lab)
            data.append(vals)
            cols.append(colors[var])
            if var not in variants:
                variants.append(var)
            x += 1.0
        centers.append((gname, (start + x - 1.0) / 2.0))
        x += 1.0
    bp = ax.boxplot(data, positions=xs, widths=0.6, patch_artist=True,
                    medianprops=dict(color="#2b2b2b", linewidth=1.6),
                    whiskerprops=dict(color="#555555"), capprops=dict(color="#555555"),
                    flierprops=dict(marker="o", markersize=4, markerfacecolor="none",
                                    markeredgecolor="#555555", alpha=0.7))
    rng = np.random.RandomState(0)
    for patch, c, xx, v in zip(bp["boxes"], cols, xs, data):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
        patch.set_edgecolor(c)
        ax.scatter(xx + rng.uniform(-0.18, 0.18, len(v)), v, s=9, color=c, alpha=0.8, zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=rotate, ha="right" if rotate else "center", fontsize=8)
    ymin = ax.get_ylim()[0]
    for g, cx in centers:
        ax.annotate(g, xy=(cx, 0), xycoords=("data", "axes fraction"), xytext=(0, -38 if not rotate else -70),
                    textcoords="offset points", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    handles = [Patch(facecolor=colors[v], alpha=0.5, edgecolor=colors[v], label=v) for v in variants]
    ax.legend(handles=handles, title=legend_title, fontsize=8, loc="best")


def save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def colormap(names):
    return {n: PALETTE[i] for i, n in enumerate(names)}


def best_so_far_matrix(runs, noncached=False):
    curves = []
    for r in runs:
        o = r["objs"][~r["cached"]] if noncached else r["objs"]
        curves.append(np.minimum.accumulate(o))
    L = max(len(c) for c in curves)
    return np.array([np.concatenate([c, np.full(L - len(c), c[-1])]) for c in curves])


def conv_panel(ax, series, colors, noncached):
    for name, runs in series:
        m = best_so_far_matrix(runs, noncached)
        x = np.arange(1, m.shape[1] + 1)
        ax.plot(x, np.median(m, axis=0), color=colors[name], label=name, linewidth=1.8)
        ax.fill_between(x, np.percentile(m, 25, axis=0), np.percentile(m, 75, axis=0),
                        color=colors[name], alpha=0.2)
    ax.set_xlabel("Evaluaciones no cacheadas (solves FMO reales)" if noncached
                  else "Número de evaluación (incluye filas cacheadas)")
    ax.set_ylabel("Mejor objetivo hasta el momento (FMO)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)


# ------------------------------------------------------------------- main
def main():
    notes = []
    N = notes.append
    N("# Análisis de experimentos nuevos (modelo FMO por órgano, K=4, CERR_Prostate)\n")
    N("Generado por `scripts/analyze_new_experiments.py`. Todas las cifras de este documento se calcularon "
      "a partir de `trajectory.csv` y `run.log` de cada corrida (15 semillas por condición). "
      "No se leyó ningún `plans_dosimetry.csv`.\n")
    N("Convenciones:\n"
      "- **Resultado final** = `FMO objective f*` de `run.log` (plan final del framework). "
      "**Mejor anytime** = mínimo de `trajectory.csv`. En búsqueda local y SA ambos coinciden "
      "(se verifica abajo); en ILS pueden diferir.\n"
      "- `trajectory.csv` **no tiene columna de tiempo** (columnas: eval, angles_deg, objective, cached, "
      "términos por estructura). La unidad de todo análisis de 'cuándo' es la **evaluación** (fila), no minutos. "
      "Las filas cacheadas no cuestan un solve de FMO, por lo que el índice de evaluación no es proporcional "
      "al tiempo de reloj; por eso también se reporta el índice contando solo evaluaciones no cacheadas.\n"
      "- Tests estadísticos: Mann-Whitney / Kruskal-Wallis / Wilcoxon bilaterales con n=15 por celda (poca potencia). "
      "Se realizan muchos tests y **no se aplicó corrección por comparaciones múltiples**; los p-valores "
      "son orientativos.\n")

    # sanity checks ------------------------------------------------------
    mism = {"ls": 0, "sa": 0, "ils": 0}
    tot = {"ls": 0, "sa": 0, "ils": 0}
    nolog = 0
    for step in STEPS:
        for cat in CATS:
            for o in ORDERS:
                for st in ("first", "best"):
                    for r in ls_runs(step, cat, o, st):
                        tot["ls"] += 1
                        mism["ls"] += abs(r["final"] - r["best"]) > 0.01
                        nolog += not r["has_final"]
            for c in CALS:
                for r in sa_runs(step, cat, c):
                    tot["sa"] += 1
                    mism["sa"] += abs(r["final"] - r["best"]) > 0.01
                    nolog += not r["has_final"]
            for p in PERTS:
                for o in ORDERS:
                    for st in ("first", "best"):
                        if st == "best" and o != "fixed":
                            continue
                        for r in ils_runs(step, cat, p, o, st):
                            tot["ils"] += 1
                            mism["ils"] += abs(r["final"] - r["best"]) > 0.01
                            nolog += not r["has_final"]
    N("## 0. Verificaciones\n")
    N(md_table(["Familia", "Corridas leídas", "Final (f*) != mínimo trayectoria (tol. 0.01)"],
               [["Búsqueda local (first+best, 3 órdenes)", tot["ls"], mism["ls"]],
                ["SA", tot["sa"], mism["sa"]],
                ["ILS (first 3 órdenes + best fijo)", tot["ils"], mism["ils"]]]))
    N("Corridas sin línea `FMO objective f*` en run.log (se usó el mínimo de la trayectoria): %d.\n" % nolog)

    # ================= A) LOCAL SEARCH =====================================
    lsdir = os.path.join(NEW, "local_search", "comparisons")
    N("## A. Búsqueda local\n")
    colors_o = colormap([ORDER_ES[o] for o in ORDERS])
    # A1
    for step in STEPS:
        fig, ax = plt.subplots(figsize=(8, 6))
        groups = [(CAT_ES[c], [(ORDER_ES[o], finals(ls_runs(step, c, o, "first")), ORDER_ES[o])
                               for o in ORDERS]) for c in CATS]
        grouped_box(ax, groups, colors_o, title="Búsqueda local, First improvement, paso %d° "
                    "(15 semillas)" % step, legend_title="Orden de escaneo")
        fig.subplots_adjust(bottom=0.18)
        save(fig, os.path.join(lsdir, "A1_first_orders_nangshift%d.png" % step))
    N("**A1.** Figuras: `local_search/comparisons/A1_first_orders_nangshift{5,10}.png` "
      "(solo First, órdenes fijo/circular/random, agrupado por catálogo).\n")

    # A2/A4 selection + tables
    chosen = {}
    a4_rows, a4_csv = [], []
    for step in STEPS:
        fig, ax = plt.subplots(figsize=(7, 6))
        groups = []
        for c in CATS:
            meds = {o: np.median(finals(ls_runs(step, c, o, "first"))) for o in ORDERS}
            co = min(meds, key=lambda k: meds[k])
            chosen[(step, c)] = co
            fr = ls_runs(step, c, co, "first")
            br = ls_runs(step, c, "fixed", "best")
            groups.append((CAT_ES[c], [("First (%s)" % ORDER_ES[co], finals(fr), "First"),
                                       ("Best (fijo)", finals(br), "Best")]))
            kwp = kw([finals(ls_runs(step, c, o, "first")) for o in ORDERS])
            mwp = mwu(finals(fr), finals(br))
            nn_idx = lambda rs: np.median([r["idx"] for r in rs])
            nnn_idx = lambda rs: np.median([r["nn_idx"] for r in rs])
            row = [step, CAT_ES[c]]
            row += [f0(meds[o]) for o in ORDERS]
            row += [fp(kwp), ORDER_ES[co], f0(np.median(finals(br))), fp(mwp),
                    f0(np.median(anytime(fr))), f0(np.median(anytime(br))),
                    "%.1f / %.1f" % (nn_idx(fr), nn_idx(br)),
                    "%.1f / %.1f" % (nnn_idx(fr), nnn_idx(br)),
                    "%.1f / %.1f" % (np.median([r["n"] for r in fr]), np.median([r["n"] for r in br]))]
            a4_rows.append(row)
            # extra: the medians of each order evals-to-final
            for o in ORDERS:
                rs = ls_runs(step, c, o, "first")
                a4_csv.append([step, c, "first", o, f0(np.median(finals(rs))),
                               np.median([r["idx"] for r in rs]), np.median([r["nn_idx"] for r in rs]),
                               np.median([r["n"] for r in rs])])
            rs = ls_runs(step, c, "fixed", "best")
            a4_csv.append([step, c, "best", "fixed", f0(np.median(finals(rs))),
                           np.median([r["idx"] for r in rs]), np.median([r["nn_idx"] for r in rs]),
                           np.median([r["n"] for r in rs])])
        grouped_box(ax, groups, colormap(["First", "Best"]), title="Búsqueda local, First vs Best, "
                    "paso %d°" % step, legend_title="Estrategia")
        fig.subplots_adjust(bottom=0.16)
        save(fig, os.path.join(lsdir, "A2_first_vs_best_nangshift%d.png" % step))
    hdr = ["Paso", "Catálogo"] + ["Mediana First-%s" % ORDER_ES[o] for o in ORDERS] + \
          ["KW p (3 órdenes)", "Orden elegido (menor mediana)", "Mediana Best (fijo)",
           "MWU p (First elegido vs Best)", "Mediana anytime First", "Mediana anytime Best",
           "Mediana eval. hasta valor final (First / Best)",
           "Mediana eval. no cacheadas hasta valor final (First / Best)",
           "Mediana evaluaciones totales (First / Best)"]
    write_csv(os.path.join(lsdir, "A4_stats.csv"), hdr, a4_rows)
    write_csv(os.path.join(lsdir, "A4_evals_to_final_per_condition.csv"),
              ["step", "catalog", "strategy", "order", "median_final", "median_eval_idx_of_final",
               "median_noncached_idx_of_final", "median_total_evals"], a4_csv)
    N("**A2.** Figuras: `A2_first_vs_best_nangshift{5,10}.png`. En cada (paso, catálogo) se compara Best (orden fijo) "
      "contra el orden de First con **menor mediana de objetivo final** (elegido con los mismos datos; "
      "sesgo de selección leve a favor de First).\n")
    N("**A4.** Estadísticos (valores = objetivo final f*; 'eval. hasta valor final' = índice de la primera fila "
      "que alcanza el mínimo de la trayectoria):\n")
    N(md_table(hdr, a4_rows))
    N("Detalle por orden en `local_search/comparisons/A4_evals_to_final_per_condition.csv`.\n")

    # A3 convergence
    for noncached in (False, True):
        suffix = "_noncached" if noncached else ""
        for step in STEPS:
            fig, axes = plt.subplots(1, 2, figsize=(13, 5))
            for ax, c in zip(axes, CATS):
                series = [(ORDER_ES[o], ls_runs(step, c, o, "first")) for o in ORDERS]
                conv_panel(ax, series, colors_o, noncached)
                ax.set_title("%s, First, paso %d° (mediana e IQR)" % (CAT_ES[c], step))
            save(fig, os.path.join(lsdir, "A3_conv_orders_nangshift%d%s.png" % (step, suffix)))
            fig, axes = plt.subplots(1, 2, figsize=(13, 5))
            for ax, c in zip(axes, CATS):
                co = chosen[(step, c)]
                series = [("First (%s)" % ORDER_ES[co], ls_runs(step, c, co, "first")),
                          ("Best (fijo)", ls_runs(step, c, "fixed", "best"))]
                conv_panel(ax, series, colormap([s[0] for s in series]), noncached)
                ax.set_title("%s, First vs Best, paso %d° (mediana e IQR)" % (CAT_ES[c], step))
            save(fig, os.path.join(lsdir, "A3_conv_first_vs_best_nangshift%d%s.png" % (step, suffix)))
    N("**A3.** Figuras de convergencia (mejor objetivo acumulado vs. índice de evaluación; línea = mediana de las "
      "15 semillas, banda = IQR; las corridas más cortas se rellenan con su valor final): "
      "`A3_conv_orders_nangshift{5,10}.png`, `A3_conv_first_vs_best_nangshift{5,10}.png` "
      "(x = número de evaluación, **incluye filas cacheadas**) y las mismas con sufijo `_noncached` "
      "(x = evaluaciones no cacheadas, es decir solves reales). First y Best requieren un número distinto de "
      "evaluaciones por iteración, por eso la comparación en x es de costo, no de iteraciones.\n")

    # ================= B) BUDGET ADEQUACY ==================================
    sadir = os.path.join(NEW, "sa", "comparisons")
    ilsdir = os.path.join(NEW, "ils", "comparisons")
    per_run = []  # family, variant, step, cat, strat, order, seed, n, nc, nn, idx, nn_idx, final, best
    for step in STEPS:
        for c in CATS:
            for cal in CALS:
                for s, r in zip(SEEDS, sa_runs(step, c, cal)):
                    per_run.append(["SA", cal, step, c, "-", "-", s, r["n"], r["nc"], r["nn"], r["idx"],
                                    r["nn_idx"], r["final"], r["best"]])
            for p in PERTS:
                for o in ORDERS:
                    for st in ("first", "best"):
                        if st == "best" and o != "fixed":
                            continue
                        for s, r in zip(SEEDS, ils_runs(step, c, p, o, st)):
                            per_run.append(["ILS", p, step, c, st, o, s, r["n"], r["nc"], r["nn"], r["idx"],
                                            r["nn_idx"], r["final"], r["best"]])
    PH = ["family", "variant", "step", "catalog", "strategy", "order", "seed", "total_evals", "cached_evals",
          "noncached_evals", "eval_idx_of_best", "noncached_idx_of_best", "final_fstar", "anytime_best"]
    write_csv(os.path.join(ilsdir, "per_run_metrics.csv"), PH, per_run)
    write_csv(os.path.join(sadir, "per_run_metrics.csv"), PH, [r for r in per_run if r[0] == "SA"])

    def sel(fam, variant=None, step=None, strat=None, order=None, cat=None):
        out = []
        for r in per_run:
            if r[0] != fam:
                continue
            if variant is not None and r[1] != variant:
                continue
            if step is not None and r[2] != step:
                continue
            if cat is not None and r[3] != cat:
                continue
            if strat is not None and r[4] != strat:
                continue
            if order is not None and r[5] != order:
                continue
            out.append(r)
        return out

    N("## B. Adecuación del presupuesto (18 min = 1080 s de reloj)\n")
    N("Para cada corrida: evaluaciones totales (filas), no cacheadas, e índice de evaluación en que se alcanzó por "
      "primera vez el mínimo de la trayectoria, como fracción del total. Como `trajectory.csv` no tiene "
      "tiempos, la fracción es **en evaluaciones**; la fracción en evaluaciones no cacheadas es una mejor "
      "aproximación del tiempo de reloj (los solves de FMO dominan el costo). "
      "Tabla agregada: SA por calendario (120 corridas = 2 pasos x 2 catálogos x 15 semillas x 2 calendarios, "
      "60 por calendario); ILS por perturbación y estrategia, agrupando pasos, catálogos y órdenes (First: 180 corridas "
      "por perturbación, Best fijo: 60).\n")
    brows, brow_csv, groups_B = [], [], []
    def bline(label, rs):
        frac = np.array([r[10] / r[7] for r in rs])
        fracn = np.array([r[11] / r[9] for r in rs])
        d, dn = dist(frac), dist(fracn)
        return [label, len(rs), f0(np.median([r[7] for r in rs])),
                f0(np.median([r[9] for r in rs])),
                "%.2f / %.2f / %.2f" % (d["med"], d["q1"], d["q3"]), "%.2f" % d["p90"],
                "%.2f / %.2f / %.2f" % (dn["med"], dn["q1"], dn["q3"]), "%.2f" % dn["p90"],
                "%.1f%%" % (100.0 * np.mean(frac <= 0.5)), "%.1f%%" % (100.0 * np.mean(fracn <= 0.5))], frac, fracn
    for cal in CALS:
        row, fr, frn = bline(CAL_ES[cal], sel("SA", cal))
        brows.append(row)
        groups_B.append((CAL_ES[cal], fr, frn))
        for step in STEPS:
            for c in CATS:
                row, _, _ = bline("%s, paso %d, %s" % (CAL_ES[cal], step, c), sel("SA", cal, step, cat=c))
                brow_csv.append(row)
    for p in PERTS:
        for st in ("first", "best"):
            lab = "ILS %s, %s" % (PERT_ES[p], st)
            row, fr, frn = bline(lab, sel("ILS", p, strat=st))
            brows.append(row)
            groups_B.append((lab, fr, frn))
            for step in STEPS:
                for c in CATS:
                    row, _, _ = bline("%s, paso %d, %s" % (lab, step, c), sel("ILS", p, step, st, cat=c))
                    brow_csv.append(row)
    BH = ["Variante", "n corridas", "Mediana evals totales", "Mediana evals no cacheadas",
          "Fracción (evals) mediana / Q1 / Q3", "Fracción p90",
          "Fracción (no cacheadas) mediana / Q1 / Q3", "Fracción (no cach.) p90",
          "% corridas con mejor en <=50% (evals)", "% corridas con mejor en <=50% (no cach.)"]
    N(md_table(BH, brows))
    write_csv(os.path.join(ilsdir, "B_budget_adequacy.csv"), BH, brows + brow_csv)
    write_csv(os.path.join(sadir, "B_budget_adequacy.csv"), BH, [r for r in brows + brow_csv if r[0].startswith("SA")])
    N("Desglose por paso y catálogo en `ils/comparisons/B_budget_adequacy.csv` y `sa/comparisons/B_budget_adequacy.csv`; "
      "métricas por corrida en `per_run_metrics.csv` de cada carpeta.\n")
    # figures (one for SA, one for ILS, ECDF both)
    for fam, pref, outdir in (("SA", "SA", sadir), ("ILS", "ILS", ilsdir)):
        gs = [g for g in groups_B if g[0].startswith(pref)]
        fig, axes = plt.subplots(1, 2, figsize=(max(10, 1.0 * len(gs) + 6), 5.5), sharey=True)
        for ax, idx, ttl in ((axes[0], 1, "Fracción del total de evaluaciones (incluye cacheadas)"),
                             (axes[1], 2, "Fracción de evaluaciones no cacheadas")):
            data = [g[idx] for g in gs]
            bp = ax.boxplot(data, patch_artist=True, medianprops=dict(color="#2b2b2b", linewidth=1.5))
            for i, (patch, g) in enumerate(zip(bp["boxes"], gs)):
                patch.set_facecolor(PALETTE[i % 10])
                patch.set_alpha(0.4)
            ax.set_xticklabels([g[0].replace("ILS ", "") for g in gs], rotation=45, ha="right", fontsize=8)
            ax.axhline(0.5, color="gray", linestyle="--", linewidth=1)
            ax.set_title(ttl, fontsize=10)
            ax.grid(axis="y", alpha=0.3)
        axes[0].set_ylabel("Posición de la evaluación donde se alcanza el mejor / total")
        fig.suptitle("%s: ¿en qué fracción del presupuesto se alcanza el mejor objetivo?" % fam)
        save(fig, os.path.join(outdir, "B_budget_fraction_box.png"))
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        for ax, idx in ((axes[0], 1), (axes[1], 2)):
            for i, g in enumerate(gs):
                v = np.sort(g[idx])
                ax.step(v, np.arange(1, len(v) + 1) / float(len(v)), where="post", color=PALETTE[i % 10],
                        label=g[0].replace("ILS ", ""))
            ax.set_xlabel("Fracción del presupuesto (evals)" if idx == 1 else "Fracción (evals no cacheadas)")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("ECDF")
        axes[1].legend(fontsize=7)
        fig.suptitle("%s: ECDF de la fracción en que se alcanza el mejor" % fam)
        save(fig, os.path.join(outdir, "B_budget_fraction_ecdf.png"))
    N("Figuras: `sa/comparisons/B_budget_fraction_{box,ecdf}.png`, `ils/comparisons/B_budget_fraction_{box,ecdf}.png`.\n")

    # ================= C) CONSISTENCY ======================================
    N("## C. Consistencia de evaluaciones con el mismo presupuesto de 18 min\n")
    N("Mediana [mín – máx] por corrida, agrupando catálogos (y órdenes en ILS First). "
      "'Cacheadas' = filas con `cached=true` (reevaluaciones de un plan ya visto). "
      "'No cacheadas' = solves FMO reales.\n")
    crows, ccsv = [], []
    def cfmt(rs, i):
        a = np.array([r[i] for r in rs])
        return "%.0f [%d – %d]" % (np.median(a), a.min(), a.max())
    def cline(lab, step, rs):
        tt = np.array([r[7] for r in rs]).astype(float)
        cc = np.array([r[8] for r in rs]).astype(float)
        return [lab, step, len(rs), cfmt(rs, 7), cfmt(rs, 8), cfmt(rs, 9), "%.1f%%" % (100 * np.median(cc / tt))]
    for step in STEPS:
        for cal in CALS:
            crows.append(cline(CAL_ES[cal], step, sel("SA", cal, step)))
        for p in PERTS:
            for st in ("first", "best"):
                crows.append(cline("ILS %s, %s" % (PERT_ES[p], st), step, sel("ILS", p, step, st)))
    CH = ["Familia / variante", "Paso", "n", "Total evals", "Cacheadas", "No cacheadas", "Mediana % cacheadas"]
    N(md_table(CH, crows))
    write_csv(os.path.join(ilsdir, "C_eval_counts.csv"), CH, crows)
    write_csv(os.path.join(sadir, "C_eval_counts.csv"), CH, [r for r in crows if r[0].startswith("SA")])
    # drivers of ILS excess: by order, and stepwise
    N("Efecto del orden de escaneo en ILS First (mediana de evals totales / no cacheadas, agrupando perturbaciones, "
      "catálogos y ambos pasos):\n")
    orows = []
    for o in ORDERS:
        rs = sel("ILS", strat="first", order=o)
        orows.append([ORDER_ES[o], len(rs), f0(np.median([r[7] for r in rs])), f0(np.median([r[9] for r in rs]))])
    N(md_table(["Orden", "n", "Mediana evals totales", "Mediana no cacheadas"], orows))
    # LS reference: evals of local search (converges, no time budget)
    lrows = []
    for step in STEPS:
        for st in ("first", "best"):
            rs = [r for c in CATS for r in ls_runs(step, c, "fixed", st)]
            lrows.append(["Búsqueda local " + st + " (fijo)", step, len(rs), f0(np.median([r["n"] for r in rs])),
                          f0(np.median([r["nc"] for r in rs])), f0(np.median([r["nn"] for r in rs]))])
    N("Referencia: búsqueda local (termina por convergencia, no por presupuesto):\n")
    N(md_table(["Familia", "Paso", "n", "Mediana total", "Mediana cacheadas", "Mediana no cacheadas"], lrows))
    sa_nn = np.median([r[9] for r in sel("SA")])
    ils_nn = np.median([r[9] for r in sel("ILS")])
    sa_t = np.median([r[7] for r in sel("SA")])
    ils_t = np.median([r[7] for r in sel("ILS")])
    N("Resumen global: mediana SA = %.0f evals totales (%.0f no cacheadas); mediana ILS = %.0f totales (%.0f no cacheadas); "
      "razón de medianas ILS/SA: %.2f (totales), %.2f (no cacheadas).\n" % (sa_t, sa_nn, ils_t, ils_nn,
                                                                           ils_t / sa_t, ils_nn / sa_nn))

    # ================= D) ILS vs SA ========================================
    N("## D. ILS vs SA con el mismo presupuesto\n")
    N("**Elección del mejor ILS:** para cada (paso, catálogo) se toma, entre las 20 configuraciones ILS "
      "(5 perturbaciones x [First con 3 órdenes + Best fijo]), la de **menor mediana de objetivo final f***. "
      "La elección usa los mismos datos de la comparación, por lo que el 'mejor ILS' está sesgado a favor "
      "(selección entre 20 configuraciones); los p-valores contra SA son optimistas para ILS. "
      "Referencia adicional: búsqueda local First con el orden elegido en A.\n")
    drows, dcsv = [], []
    for step in STEPS:
        fig, ax = plt.subplots(figsize=(9, 6.5))
        groups = []
        for c in CATS:
            cand = {}
            for p in PERTS:
                for o in ORDERS:
                    cand[(p, o, "first")] = finals(ils_runs(step, c, p, o, "first"))
                cand[(p, "fixed", "best")] = finals(ils_runs(step, c, p, "fixed", "best"))
            bk = min(cand, key=lambda k: np.median(cand[k]))
            ils_v = cand[bk]
            sl, sg = finals(sa_runs(step, c, "linear")), finals(sa_runs(step, c, "geometric"))
            lsv = finals(ls_runs(step, c, chosen[(step, c)], "first"))
            lab = "ILS %s/%s/%s" % (bk[0], ORDER_ES[bk[1]], bk[2])
            default = cand[("range3", "fixed", "first")]
            groups.append((CAT_ES[c], [("ILS mejor\n(%s, %s, %s)" % (bk[0], ORDER_ES[bk[1]], bk[2]), ils_v, "ILS mejor"),
                                       ("SA lineal", sl, "SA lineal"), ("SA geom.", sg, "SA geométrico"),
                                       ("LS First\n(%s)" % ORDER_ES[chosen[(step, c)]], lsv, "LS First (ref.)")]))
            drows.append([step, CAT_ES[c], lab, f0(np.median(ils_v)), f0(np.median(sl)), f0(np.median(sg)),
                          f0(np.median(lsv)), f0(np.median(default)),
                          fp(mwu(ils_v, sl)), fp(mwu(ils_v, sg)), fp(mwu(sl, sg)), fp(mwu(ils_v, lsv)),
                          fp(mwu(default, sl)), fp(mwu(default, sg))])
        grouped_box(ax, groups, colormap(["ILS mejor", "SA lineal", "SA geométrico", "LS First (ref.)"]),
                    title="ILS (mejor config.) vs SA vs búsqueda local, paso %d°" % step,
                    legend_title="Algoritmo")
        fig.subplots_adjust(bottom=0.2)
        save(fig, os.path.join(ilsdir, "ils_vs_sa_nangshift%d.png" % step))
    DH = ["Paso", "Catálogo", "Mejor ILS (pert/orden/estrategia)", "Mediana ILS mejor", "Mediana SA lineal",
          "Mediana SA geom.", "Mediana LS First", "Mediana ILS range3-First-fijo",
          "MWU p ILS mejor vs SA lin", "MWU p ILS mejor vs SA geom", "MWU p SA lin vs geom",
          "MWU p ILS mejor vs LS First", "MWU p ILS range3-fijo vs SA lin", "MWU p ILS range3-fijo vs SA geom"]
    N(md_table(DH, drows))
    write_csv(os.path.join(ilsdir, "D_ils_vs_sa.csv"), DH, drows)
    N("Figuras: `ils/comparisons/ils_vs_sa_nangshift{5,10}.png`. Se añade la columna con ILS range3 (valor de "
      "literatura), First orden fijo, sin selección, como comparación libre de sesgo de selección.\n")

    # ================= E) ILS perturbation ================================
    N("## E. Tamaño y tipo de perturbación en ILS\n")
    pert_colors = colormap([PERT_ES[p] for p in PERTS])
    erows, pcsv = [], []
    for step in STEPS:
        for st, orders in (("first", ["fixed"]), ("best", ["fixed"])):
            fig, ax = plt.subplots(figsize=(12, 6))
            groups = [(CAT_ES[c], [(PERT_ES[p], finals(ils_runs(step, c, p, "fixed", st)), PERT_ES[p])
                                   for p in PERTS]) for c in CATS]
            grouped_box(ax, groups, pert_colors, title="ILS %s (orden fijo), paso %d°: perturbaciones" % (st, step),
                        legend_title="Perturbación", rotate=30)
            fig.subplots_adjust(bottom=0.28)
            save(fig, os.path.join(ilsdir, "E_perturbation_%s_nangshift%d.png" % (st, step)))
        # order effect within perturbation (first)
        fig, axes = plt.subplots(2, 1, figsize=(12, 9))
        for ax, c in zip(axes, CATS):
            groups = [(PERT_ES[p], [(ORDER_ES[o], finals(ils_runs(step, c, p, o, "first")), ORDER_ES[o])
                                    for o in ORDERS]) for p in PERTS]
            grouped_box(ax, groups, colors_o, title="ILS First, %s, paso %d°: orden de escaneo por perturbación"
                        % (CAT_ES[c], step), legend_title="Orden")
        fig.subplots_adjust(hspace=0.45)
        save(fig, os.path.join(ilsdir, "E_order_by_perturbation_first_nangshift%d.png" % step))
    # tables
    kwrows, mwrows, ordrows = [], [], []
    for step in STEPS:
        for c in CATS:
            for st in ("first", "best"):
                for o in (ORDERS if st == "first" else ["fixed"]):
                    v = {p: finals(ils_runs(step, c, p, o, st)) for p in PERTS}
                    kwrows.append([step, CAT_ES[c], st, ORDER_ES[o]] + [f0(np.median(v[p])) for p in PERTS] +
                                  [fp(kw([v[p] for p in PERTS]))])
                    mwrows.append([step, CAT_ES[c], st, ORDER_ES[o],
                                   fp(mwu(v["range3"], v["mag2"])), fp(mwu(v["mag2"], v["mag2-stepwise"])),
                                   fp(mwu(v["adaptive"], v["adaptivemag-stepwise"])),
                                   fp(mwu(v["range3"], v["adaptive"]))])
            for p in PERTS:
                v = {o: finals(ils_runs(step, c, p, o, "first")) for o in ORDERS}
                ordrows.append([step, CAT_ES[c], PERT_ES[p]] + [f0(np.median(v[o])) for o in ORDERS] +
                               [fp(kw([v[o] for o in ORDERS]))])
    EH = ["Paso", "Catálogo", "Estrategia", "Orden"] + ["Mediana " + PERT_ES[p] for p in PERTS] + ["KW p (5 pert.)"]
    N("**E1.** Medianas del objetivo final por perturbación y Kruskal-Wallis entre las 5 perturbaciones:\n")
    N(md_table(EH, kwrows))
    write_csv(os.path.join(ilsdir, "E1_perturbation_kw.csv"), EH, kwrows)
    MH = ["Paso", "Catálogo", "Estrategia", "Orden", "MWU p range3 vs mag2", "MWU p mag2 vs mag2-stepwise",
          "MWU p adaptive vs adaptivemag-stepwise", "MWU p range3 vs adaptive"]
    N("**E2.** Comparaciones entre perturbaciones (Mann-Whitney no pareado). 'stepwise vs no stepwise' se mide con "
      "mag2 vs mag2-stepwise (misma magnitud, 2 ángulos) y adaptive vs adaptivemag-stepwise (esta última cambia también "
      "el tipo de magnitud, por lo que no es un contraste limpio de stepwise):\n")
    N(md_table(MH, mwrows))
    write_csv(os.path.join(ilsdir, "E2_perturbation_pairwise.csv"), MH, mwrows)
    OH = ["Paso", "Catálogo", "Perturbación"] + ["Mediana " + ORDER_ES[o] for o in ORDERS] + ["KW p (3 órdenes)"]
    same = tot_s = 0
    for step in STEPS:
        for c in CATS:
            for st in ("first", "best"):
                for o in (ORDERS if st == "first" else ["fixed"]):
                    a, b = finals(ils_runs(step, c, "mag2", o, st)), finals(ils_runs(step, c, "mag2-stepwise", o, st))
                    same += int((np.abs(a - b) < 0.01).sum())
                    tot_s += len(a)
    N("Semillas en que mag2 y mag2-stepwise terminan con el mismo f* (|dif|<0.01): %d de %d (%.1f%%); "
      "las medianas y p-valores casi idénticos entre ambas columnas reflejan esto.\n" % (same, tot_s, 100.0 * same / tot_s))
    N("**E3.** Efecto del orden de escaneo dentro de cada perturbación (ILS First):\n")
    N(md_table(OH, ordrows))
    write_csv(os.path.join(ilsdir, "E3_order_within_perturbation.csv"), OH, ordrows)
    N("Figuras: `ils/comparisons/E_perturbation_{first,best}_nangshift{5,10}.png`, "
      "`E_order_by_perturbation_first_nangshift{5,10}.png`.\n")
    n_sig = sum(1 for r in kwrows if r[-1] != "n/d" and (r[-1].startswith("<") or float(r[-1]) < 0.05))
    N("De %d tests KW entre perturbaciones, %d tienen p<0.05 (sin corrección por comparaciones múltiples).\n"
      % (len(kwrows), n_sig))

    # pairing check
    N("**E4. Pareo ILS vs búsqueda local por semilla.**\n")
    ini = {}
    for fam in ("ls", "ils"):
        eqc = totc = 0
        for step in STEPS:
            for c in CATS:
                for sd in range(15):
                    if fam == "ls":
                        vals = [ls_runs(step, c, o, st)[sd]["first_obj"] for o in ORDERS for st in ("first", "best")]
                    else:
                        vals = [ils_runs(step, c, p, o, "first")[sd]["first_obj"] for p in PERTS for o in ORDERS]
                    totc += 1
                    eqc += int(max(vals) - min(vals) < 1e-6)
        ini[fam] = (eqc, totc)
    N("Control: la solución inicial (fila 1) es la misma entre todas las configuraciones de una misma familia para "
      "cada (paso, catálogo, semilla): búsqueda local %d/%d, ILS %d/%d. Es decir, la semilla fija el inicio dentro de "
      "cada familia, pero no entre familias.\n" % (ini["ls"][0], ini["ls"][1], ini["ils"][0], ini["ils"][1]))
    eq_all, tot_all = 0, 0
    pair_rows, pair_csv = [], []
    for step in STEPS:
        for c in CATS:
            for st in ("first", "best"):
                for o in (ORDERS if st == "first" else ["fixed"]):
                    ls_r = ls_runs(step, c, o, st)
                    for p in PERTS:
                        il = ils_runs(step, c, p, o, st)
                        eq = sum(abs(a["first_obj"] - b["first_obj"]) < 1e-6 for a, b in zip(il, ls_r))
                        eq_all += eq
                        tot_all += len(il)
                        valid = eq == len(il)
                        a, b = finals(il), finals(ls_r)
                        d = a - b
                        if valid:
                            pv, test = wilcoxon(d), "Wilcoxon"
                        else:
                            pv, test = mwu(a, b), "MWU (no pareado)"
                        row = [step, CAT_ES[c], st, ORDER_ES[o], PERT_ES[p], "%d/15" % eq, f0(np.median(b)),
                               f0(np.median(a)), f0(np.median(d)), "%d/15" % int((d < 0).sum()), test, fp(pv)]
                        pair_csv.append(row)
                        if o == "fixed":
                            pair_rows.append(row)
    N("Comparación del objetivo de la fila 1 de `trajectory.csv` (solución inicial) entre ILS y búsqueda local "
      "para igual (paso, catálogo, semilla, estrategia, orden): coincidencias %d de %d (%.1f%%).\n"
      % (eq_all, tot_all, 100.0 * eq_all / tot_all))
    PHD = ["Paso", "Catálogo", "Estrategia", "Orden", "Perturbación", "Semillas con igual inicio",
           "Mediana LS final", "Mediana ILS final", "Mediana (ILS - LS)", "Semillas con ILS<LS", "Test", "p"]
    if eq_all == tot_all:
        N("El pareo es válido (mismo inicio en todas las semillas). Mejora pareada = ILS final - LS final "
          "(negativo = ILS mejor), test de Wilcoxon de rangos con signo. Se muestran las filas de orden fijo; "
          "todas en `E4_ils_vs_ls_paired.csv`.\n")
    else:
        N("El pareo NO es válido en todas las comparaciones; se usa Mann-Whitney no pareado donde no coincide el "
          "inicio (columna 'Test'). Filas de orden fijo; todas en `E4_ils_vs_ls_paired.csv`.\n")
    N(md_table(PHD, pair_rows))
    write_csv(os.path.join(ilsdir, "E4_ils_vs_ls_paired.csv"), PHD, pair_csv)

    # extra: framework final vs anytime best (ILS)
    ils_gap = [(r[12], r[13]) for r in sel("ILS")]
    nworse = sum(1 for a, b in ils_gap if a - b > 0.01)
    N("## Anexo: resultado final f* vs mejor anytime en ILS\n")
    N("Corridas ILS donde f* es peor que el mínimo de la trayectoria (>0.01): %d de %d. "
      "Mediana de la diferencia entre ellas: %s.\n" % (
          nworse, len(ils_gap),
          f1(np.median([a - b for a, b in ils_gap if a - b > 0.01])) if nworse else "n/a"))

    arows = []
    for st in ("first", "best"):
        for p in PERTS:
            rs = sel("ILS", p, strat=st)
            w = [r[12] - r[13] for r in rs if r[12] - r[13] > 0.01]
            arows.append([st, PERT_ES[p], len(rs), len(w), f1(np.median(w)) if w else "n/a"])
    N("Por estrategia y perturbación:\n")
    N(md_table(["Estrategia", "Perturbación", "n", "f* peor que anytime", "Mediana de la diferencia"], arows))
    path = os.path.join(NEW, "ANALYSIS_NOTES.md")
    with open(path, "w") as f:
        f.write("\n".join(notes))
    print("wrote", path)


if __name__ == "__main__":
    main()
