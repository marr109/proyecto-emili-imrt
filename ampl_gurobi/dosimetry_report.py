"""Dosimetric quality of a plan: dose per structure for a given angle set.
Uses the exact FMO solution (cerr_fmo_osqp) and reports, per structure:
Dmin, Dmean, Dmax, and the clinically relevant indicators.
Limits are the PROVISIONAL ones used by the solver (dmin 65, dmax 50 Gy);
they are not a clinical prescription and must be replaced before any clinical claim."""
import sys
import numpy as np
from cerr_fmo_ampl import CerrInstance, PTV_ORGANS, OAR_ORGANS, DMIN, DMAX_OAR
from solvers import build_and_solve

def report(angles, w_over_organ, inst=None):
    inst = inst or CerrInstance()
    out = build_and_solve(inst, angles, w_over_organ)
    if out["status"] not in ("solved", "solved inaccurate"):
        raise RuntimeError(f"solver status {out['status']}")
    D_ptv, D_oar = out["D_ptv"], out["D_oar"]
    ptv_split = inst.n_rows[PTV_ORGANS[0]]
    organ_rows = {"PTVHD": D_ptv[:ptv_split], "PTVLD": D_ptv[ptv_split:]}
    b_rows = inst.n_rows[OAR_ORGANS[0]]
    organ_rows["BLADDER"] = D_oar[:b_rows]
    organ_rows["RECTUM"] = D_oar[b_rows:]
    rows = []
    for name, d in organ_rows.items():
        d = np.asarray(d)
        is_ptv = name.startswith("PTV")
        rows.append({
            "structure": name,
            "voxels": len(d),
            "Dmin": float(d.min()),
            "Dmean": float(d.mean()),
            "Dmax": float(d.max()),
            "D95": float(np.percentile(d, 5)) if is_ptv else None,   # dose covering 95% of volume
            "D2": float(np.percentile(d, 98)) if is_ptv else None,   # near-max (2% of volume)
            "V_limit_pct": float(100 * np.mean(d >= DMIN)) if is_ptv else float(100 * np.mean(d > DMAX_OAR)),
        })
    return {"angles": sorted(angles), "objective": out["objective"], "structures": rows}

def print_report(r):
    print(f"angles={r['angles']}  objective={r['objective']:.4f}")
    print(f"{'structure':<9}{'vox':>7}{'Dmin':>9}{'Dmean':>9}{'Dmax':>9}{'D95':>9}{'D2':>9}  note")
    for s in r["structures"]:
        f = lambda v: f"{v:9.2f}" if v is not None else f"{'-':>9}"
        if s["structure"].startswith("PTV"):
            note = f"% vol >= {DMIN} Gy: {s['V_limit_pct']:.1f}"
        else:
            note = f"% vol > {DMAX_OAR} Gy: {s['V_limit_pct']:.1f}"
        print(f"{s['structure']:<9}{s['voxels']:>7}{f(s['Dmin'])}{f(s['Dmean'])}{f(s['Dmax'])}{f(s['D95'])}{f(s['D2'])}  {note}")

if __name__ == "__main__":
    angles = [int(a) for a in sys.argv[1].split(",")] if len(sys.argv) > 1 else [35, 115, 165, 285]
    w = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.5, 0.8]
    print_report(report(angles, w))
