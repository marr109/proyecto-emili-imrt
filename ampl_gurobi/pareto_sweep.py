"""Pareto front (FMO level) for a FIXED angle set, quadratic-penalty objectives.

Objectives (mean squared slack, so they are comparable across structures):
  f_ptv     PTV underdose        f_bladder  bladder overdose      f_rectum  rectum overdose
Weighted-sum sweep over the simplex (a, b, c) = weights of (PTV, bladder, rectum).
The problem is convex, so every Pareto-optimal plan is reachable by some weighting.
Hot-spot term stays fixed (W_HOT) and is NOT one of the swept objectives.

Usage: python pareto_sweep.py 25,105,160,280 [divisions]
Env:   EMILI_DMIN (prescription, default 65), OUT_TAG
"""
import os, sys, json, time, itertools
from cerr_fmo_ampl import CerrInstance
from solvers import build_and_solve

W_HOT = 0.1
FLOOR = 0.02            # keeps a zero-weighted objective from being ignored completely

def simplex(div):
    for i in range(div + 1):
        for j in range(div + 1 - i):
            yield (i / div, j / div, (div - i - j) / div)

def dominated(p, others):
    return any(all(o[k] <= p[k] for k in range(3)) and any(o[k] < p[k] for k in range(3)) for o in others)

if __name__ == "__main__":
    angles = sorted(int(a) for a in sys.argv[1].split(","))
    div = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    tag = os.environ.get("OUT_TAG", "")
    inst = CerrInstance()
    pts = []
    for (a, b, c) in simplex(div):
        wp, wb, wr = max(a, FLOOR), max(b, FLOOR), max(c, FLOOR)
        t = time.time()
        out = build_and_solve(inst, angles, [wb, wr], w_under=wp, w_ptv_over=W_HOT, normalize=True)
        rec = {"weights": [wp, wb, wr], "status": out["status"],
               "f_ptv": out["f_ptv"], "f_bladder": out["f_bladder"], "f_rectum": out["f_rectum"]}
        if out["status"] in ("solved", "solved inaccurate"):
            import numpy as np
            nb = inst.n_rows["BLADDER"]; ps = inst.n_rows["PTVHD"]
            D_ptv, D_oar = out["D_ptv"], out["D_oar"]
            rec.update(ptv_D95=float(np.percentile(D_ptv, 5)), ptv_mean=float(D_ptv.mean()),
                       bladder_mean=float(D_oar[:nb].mean()), bladder_max=float(D_oar[:nb].max()),
                       rectum_mean=float(D_oar[nb:].mean()), rectum_max=float(D_oar[nb:].max()))
        pts.append(rec)
        print(f"w=({wp:.2f},{wb:.2f},{wr:.2f}) f=({rec['f_ptv']:.3f},{rec['f_bladder']:.3f},{rec['f_rectum']:.3f}) "
              f"{out['status']} ({time.time()-t:.0f}s)", flush=True)
    ok = [p for p in pts if p["status"] in ("solved", "solved inaccurate")]
    objs = [(p["f_ptv"], p["f_bladder"], p["f_rectum"]) for p in ok]
    for p, o in zip(ok, objs):
        p["nondominated"] = not dominated(o, [x for x in objs if x != o])
    json.dump({"angles": angles, "dmin": float(os.environ.get("EMILI_DMIN", 65)), "points": pts},
              open(f"pareto_{'-'.join(map(str, angles))}{tag}.json", "w"), indent=2)
    print("nondominated:", sum(p["nondominated"] for p in ok), "of", len(ok))
