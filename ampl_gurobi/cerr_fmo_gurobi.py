"""CERR_Prostate FMO with per-organ OAR weights, solved with AMPL + Gurobi (production path).

Same interface and return keys as cerr_fmo_osqp.build_and_solve, so the search,
Pareto sweep and dosimetry scripts are solver-agnostic (see solvers.py).
"""
import os, tempfile
import numpy as np
import scipy.sparse as sp
from amplpy import AMPL
from cerr_fmo_ampl import CerrInstance, MOD, PTV_ORGANS, OAR_ORGANS

def _dose_matrices(inst, active):
    dim_of, n = {}, 0
    for a in active:
        s, e = inst.beamlet_range[a]
        for g in range(s - 1, e):
            dim_of[g] = n; n += 1
    n_ptv = sum(inst.n_rows[o] for o in PTV_ORGANS)
    n_oar = sum(inst.n_rows[o] for o in OAR_ORGANS)
    P, O = ([], [], []), ([], [], [])
    for a in active:
        data = inst.angle(a); g0 = inst.beamlet_range[a][0] - 1
        for group, acc in ((PTV_ORGANS, P), (OAR_ORGANS, O)):
            for organ in group:
                off = inst.row_off[organ]
                for v, b, d in data[organ]:
                    acc[0].append(off + inst.voxel_row[organ][v]); acc[1].append(dim_of[g0 + b - 1]); acc[2].append(d)
    A_ptv = sp.csr_matrix((P[2], (P[0], P[1])), shape=(n_ptv, n))
    A_oar = sp.csr_matrix((O[2], (O[0], O[1])), shape=(n_oar, n))
    return A_ptv, A_oar, dim_of, n_ptv, n_oar

def build_and_solve(inst, active, w_over_organ, w_under=None, w_ptv_over=None, normalize=False, **_):
    active = sorted(active)
    dat = inst.dat_text(active, w_over_organ, w_under=w_under, w_ptv_over=w_ptv_over)
    ampl = AMPL(); ampl.option["solver"] = "gurobi"; ampl.option["show_stats"] = 0
    ampl.read(MOD)
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
        f.write(dat); name = f.name
    ampl.read_data(name); os.unlink(name)
    ampl.get_parameter("normalize").set(1 if normalize else 0)
    ampl.solve()
    res = ampl.solve_result
    if res != "solved":
        ampl.close(); return {"status": res, "objective": float("inf")}
    A_ptv, A_oar, dim_of, n_ptv, n_oar = _dose_matrices(inst, active)
    x_d = ampl.get_variable("x").to_dict()
    x = np.zeros(len(dim_of))
    for g, i in dim_of.items(): x[i] = x_d[g]
    u_d, v_d = ampl.get_variable("u").to_dict(), ampl.get_variable("v").to_dict()
    u = np.array([u_d[i] for i in range(n_ptv)]); v = np.array([v_d[i] for i in range(n_oar)])
    obj = ampl.get_objective("fmo_objective").value()
    ampl.close()
    org = np.array([0] * inst.n_rows[OAR_ORGANS[0]] + [1] * inst.n_rows[OAR_ORGANS[1]])
    return {"status": "solved", "objective": obj, "f_ptv": float(np.mean(u ** 2)),
            "f_bladder": float(np.mean(v[org == 0] ** 2)), "f_rectum": float(np.mean(v[org == 1] ** 2)),
            "x": x, "u": u, "v": v, "D_ptv": A_ptv @ x, "D_oar": A_oar @ x}
