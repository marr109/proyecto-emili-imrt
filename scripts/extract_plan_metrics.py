#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrae las metricas dosimetricas que reportPlan() imprime al cierre de cada
corrida y las emite como CSV o como tabla LaTeX.

El bloque clinico del run.log tiene tres partes: estadisticas por organo en Gy,
restricciones DVH con su estado, e indices de calidad de plan. Este script lee
las tres y las aplana a una fila por corrida.

ADVERTENCIA sobre la interpretacion: mientras w_ptv_over valga 0 el modelo no
penaliza la sobredosis en el PTV, de modo que D2, Dmax, CI y HI describen una
distribucion que nunca se pidio homogenea. No son medidas de calidad de plan
bajo esa configuracion. Ver imrt/cerr_instance.cpp.

Uso:
  python3 scripts/extract_plan_metrics.py --log a.log --log b.log --format latex
  python3 scripts/extract_plan_metrics.py --glob 'experiments/sa/**/run.log'
"""
import argparse, glob, io, os, re, sys

ORG_RE = re.compile(
    r"^(\w+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s*$")
IDX_RE = {
    "CI":  re.compile(r"Conformity Index.*?=\s*([\d.]+)"),
    "HI":  re.compile(r"Homogeneity Index.*?=\s*([\d.]+)"),
    "V95": re.compile(r"PTV coverage.*?=\s*([\d.]+)\s*%"),
    "obj": re.compile(r"FMO objective f\*\s*:\s*([\d.]+)"),
    "Rx":  re.compile(r"Prescription dose Rx\s*:\s*([\d.]+)"),
}


def parse(path):
    txt = io.open(path, encoding="utf-8", errors="replace").read()
    if "Conformity Index" not in txt:
        return None
    rec = {"log": path, "organs": {}}
    for k, rx in IDX_RE.items():
        m = rx.search(txt)
        rec[k] = float(m.group(1)) if m else None
    # Las estadisticas por organo viven entre su encabezado y el bloque DVH.
    blk = txt.split("Per-organ dose statistics")
    if len(blk) > 1:
        for line in blk[1].split("DVH constraints")[0].splitlines():
            m = ORG_RE.match(line.strip())
            if m:
                g = m.groups()
                rec["organs"][g[0]] = dict(
                    voxels=int(g[1]), dmin=float(g[2]), dmean=float(g[3]),
                    dmax=float(g[4]),
                    d95=(None if g[5] == "-" else float(g[5])),
                    d5=float(g[6]), d2=float(g[7]))
    rec["viol"] = len(re.findall(r"\bVIOL\b", txt))
    rec["ok"] = len(re.findall(r"\bOK\b", txt))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", action="append", default=[])
    ap.add_argument("--glob")
    ap.add_argument("--label", action="append", default=[],
                    help="etiqueta por log, en el mismo orden")
    ap.add_argument("--organ", default="PTVHD", help="organo a tabular")
    ap.add_argument("--format", choices=("csv", "latex", "text"), default="text")
    a = ap.parse_args()

    paths = list(a.log)
    if a.glob:
        paths += sorted(glob.glob(a.glob, recursive=True))
    if not paths:
        print("sin logs", file=sys.stderr)
        return 1

    rows = []
    for i, p in enumerate(paths):
        r = parse(p)
        if r is None:
            print("sin bloque clinico: %s" % p, file=sys.stderr)
            continue
        r["label"] = a.label[i] if i < len(a.label) else os.path.basename(p)
        rows.append(r)
    if not rows:
        return 1

    org = a.organ
    cols = ["label", "obj", "dmean", "dmax", "d95", "d2", "CI", "HI", "V95", "viol"]
    def cell(r, c):
        if c == "label": return r["label"]
        if c in ("obj", "CI", "HI", "V95", "viol"): return r.get(c)
        return r["organs"].get(org, {}).get(c)

    if a.format == "csv":
        print(",".join(cols))
        for r in rows:
            print(",".join("" if cell(r, c) is None else str(cell(r, c)) for c in cols))
    elif a.format == "latex":
        hdr = ["Configuración", "$f^*$", "Dmean", "Dmax", "D95", "D2", "CI", "HI", "V95\\,\\%", "VIOL"]
        print("\\begin{tabular}{|l|" + "c|" * (len(hdr) - 1) + "}")
        print("\\hline")
        print(" & ".join("\\textbf{%s}" % h for h in hdr) + " \\\\")
        print("\\hline")
        for r in rows:
            out = [r["label"]]
            for c in cols[1:]:
                v = cell(r, c)
                if v is None: out.append("---")
                elif c == "viol": out.append("%d" % v)
                elif c == "obj": out.append(("%.1f" % v).replace(".", "{,}"))
                else: out.append(("%.2f" % v).replace(".", "{,}"))
            print(" & ".join(out) + " \\\\")
            print("\\hline")
        print("\\end{tabular}")
    else:
        w = "%-26s %12s %9s %9s %9s %9s %7s %7s %8s %6s"
        print(w % ("configuracion", "objetivo", "Dmean", "Dmax", "D95", "D2", "CI", "HI", "V95%", "VIOL"))
        for r in rows:
            vals = []
            for c in cols[1:]:
                v = cell(r, c)
                vals.append("---" if v is None else ("%d" % v if c == "viol" else "%.2f" % v))
            print(w % tuple([r["label"][:26]] + vals))
    return 0


if __name__ == "__main__":
    sys.exit(main())
