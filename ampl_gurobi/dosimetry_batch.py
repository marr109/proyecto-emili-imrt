"""Dosimetric metrics of the final plan of every local-search run.

Usage:
    python dosimetry_batch.py <experiments_root> [--out plans_dosimetry.csv] [--geud]

Scans <experiments_root> for run.log files (layout
nangshift<step>/<catalog>[-suffix]/data/<strategy>/seedNN/run.log), reads the
final angles ("Selected angles (deg)"), re-solves the FMO exactly for each
distinct angle set and writes one CSV row per run.

Model settings come from the same environment variables as the C++ solver and
the other scripts (EMILI_DMIN, EMILI_W_PTV_OVER, EMILI_EXCLUDE_OVERLAP,
EMILI_W_OVER_ORGANS="bladder,rectum"), so they must match the batch run.
--geud adds gEUD columns (a = -10 PTV, 8 rectum, 2 bladder, the values used in
the Cabrera-Guerrero papers). They are for comparison with those papers only;
the model itself does not optimize gEUD.
"""
import argparse
import csv
import json
import os
import re
import sys

import numpy as np

from cerr_fmo_ampl import (CerrInstance, PTV_ORGANS, OAR_ORGANS, DMIN, DMAX_OAR,
                           DMAX_PTV, W_UNDER, W_PTV_OVER)
from solvers import build_and_solve

ANGLES_RE = re.compile(r"Selected angles \(deg\)\s*:\s*\[([^\]]*)\]")
GEUD_A = {"ptv": -10.0, "rectum": 8.0, "bladder": 2.0}


def geud(d, a):
    d = np.maximum(np.asarray(d, dtype=float), 1e-9)
    return float(np.mean(d ** a) ** (1.0 / a))


def find_runs(root):
    runs = []
    for dirpath, _, files in os.walk(root):
        if "run.log" not in files:
            continue
        with open(os.path.join(dirpath, "run.log"), errors="replace") as f:
            m = ANGLES_RE.search(f.read())
        if not m:
            continue
        angles = sorted(int(a) for a in m.group(1).replace(",", " ").split())
        runs.append((os.path.relpath(dirpath, root), angles))
    return sorted(runs)


def metrics(inst, angles, w_organs, with_geud):
    out = build_and_solve(inst, angles, w_organs)
    if out["status"] not in ("solved", "solved inaccurate"):
        return {"status": out["status"]}
    Dp, Do = np.asarray(out["D_ptv"]), np.asarray(out["D_oar"])
    nb = inst.n_rows[OAR_ORGANS[0]]
    bladder, rectum = Do[:nb], Do[nb:]
    u = np.maximum(DMIN - Dp, 0.0)
    w = np.maximum(Dp - DMAX_PTV, 0.0)
    row = {
        "status": out["status"],
        "objective": out["objective"],
        "term_ptv_under": W_UNDER * float(np.sum(u ** 2)),
        "term_ptv_hot": W_PTV_OVER * float(np.sum(w ** 2)),
        "term_bladder": w_organs[0] * float(np.sum(np.maximum(bladder - DMAX_OAR, 0) ** 2)),
        "term_rectum": w_organs[1] * float(np.sum(np.maximum(rectum - DMAX_OAR, 0) ** 2)),
        "ptv_dmin": float(Dp.min()), "ptv_dmean": float(Dp.mean()),
        "ptv_d95": float(np.percentile(Dp, 5)), "ptv_d2": float(np.percentile(Dp, 98)),
        "ptv_dmax": float(Dp.max()), "ptv_pct_ge_rx": float(100 * np.mean(Dp >= DMIN)),
        "bladder_dmean": float(bladder.mean()), "bladder_dmax": float(bladder.max()),
        "bladder_pct_gt_limit": float(100 * np.mean(bladder > DMAX_OAR)),
        "rectum_dmean": float(rectum.mean()), "rectum_dmax": float(rectum.max()),
        "rectum_pct_gt_limit": float(100 * np.mean(rectum > DMAX_OAR)),
    }
    if with_geud:
        row["geud_ptv"] = geud(Dp, GEUD_A["ptv"])
        row["geud_rectum"] = geud(rectum, GEUD_A["rectum"])
        row["geud_bladder"] = geud(bladder, GEUD_A["bladder"])
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--out", default="plans_dosimetry.csv")
    ap.add_argument("--geud", action="store_true")
    args = ap.parse_args()

    w_organs = [float(x) for x in os.environ.get("EMILI_W_OVER_ORGANS", "0.5,0.8").split(",")]
    runs = find_runs(args.root)
    if not runs:
        sys.exit(f"no run.log with final angles under {args.root}")
    inst = CerrInstance()
    # Plan metrics are cached on disk next to the output so a crash (e.g. an expired
    # solver license lease) can resume without re-solving finished plans. The cache
    # is only valid for the same model settings, so delete it if they change.
    cache_path = args.out + ".cache.json"
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            cache = {tuple(int(a) for a in k.split("-")): v for k, v in json.load(f).items()}
        print(f"resuming from {cache_path}: {len(cache)} plans already solved")
    rows = []
    for i, (rel, angles) in enumerate(runs, 1):
        key = tuple(angles)
        if key not in cache:
            cache[key] = metrics(inst, angles, w_organs, args.geud)
            with open(cache_path, "w") as f:
                json.dump({"-".join(map(str, k)): v for k, v in cache.items()}, f)
        rows.append({"run": rel, "angles": "-".join(map(str, angles)), **cache[key]})
        print(f"[{i}/{len(runs)}] {rel} {angles}")
    cols = list(dict.fromkeys(c for r in rows for c in r))
    with open(args.out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)
    print(f"wrote {args.out}: {len(rows)} runs, {len(cache)} distinct plans")


if __name__ == "__main__":
    main()
