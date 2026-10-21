"""CERR_Prostate FMO with per-organ OAR weights, solved with OSQP.

Same QP as fmo_organs.mod, built directly as sparse matrices from the CERR
instance (see cerr_fmo_ampl.CerrInstance). Used because the AMPL demo license
(300 vars/constraints) and the Gurobi restricted license (2000) are too small
for the CERR_Prostate QP (~16k variables). OSQP is the solver EMILI's README
already names for this problem.

Decision variables z = [x (n beamlets), u (PTV rows), v (OAR rows), w (PTV rows)]

  min  w_under * |u|^2 + sum_b w_over[organ(b)] * v_b^2 + w_ptv_over * |w|^2
  s.t. A_ptv x + u >= dmin          (PTV underdose)
       A_oar x - v <= dmax          (OAR overdose)
       A_ptv x - w <= dmax_ptv      (PTV hot spot)
       0 <= x <= M,  u, v, w >= 0
"""
import numpy as np
import scipy.sparse as sp
import osqp

from cerr_fmo_ampl import (CerrInstance, PTV_ORGANS, OAR_ORGANS, DMIN, DMAX_OAR,
                           DMAX_PTV, MAX_INTENSITY, W_UNDER, W_PTV_OVER)


def build_and_solve(inst, active, w_over_organ, eps=1e-8, max_iter=200000,
                    w_under=None, w_ptv_over=None, normalize=False):
    w_under = W_UNDER if w_under is None else w_under
    w_ptv_over = W_PTV_OVER if w_ptv_over is None else w_ptv_over
    active = sorted(active)
    K = len(active)
    # Local beamlet index for each active angle, in the order of the dimlets.
    dim_of = {}
    n = 0
    for a in active:
        s, e = inst.beamlet_range[a]
        for g in range(s - 1, e):
            dim_of[g] = n
            n += 1

    n_ptv = inst.n_rows[PTV_ORGANS[0]] + inst.n_rows[PTV_ORGANS[1]]
    n_oar = inst.n_rows[OAR_ORGANS[0]] + inst.n_rows[OAR_ORGANS[1]]

    ri, cj, vals = [], [], []
    rp, cp, vp = [], [], []          # PTV dose entries
    ro, co, vo = [], [], []          # OAR dose entries
    for a in active:
        data = inst.angle(a)
        g0 = inst.beamlet_range[a][0] - 1
        for organ in PTV_ORGANS:
            off = inst.row_off[organ]
            for v, b, d in data[organ]:
                rp.append(off + inst.voxel_row[organ][v])
                cp.append(dim_of[g0 + b - 1])
                vp.append(d)
        for organ in OAR_ORGANS:
            off = inst.row_off[organ]
            for v, b, d in data[organ]:
                ro.append(off + inst.voxel_row[organ][v])
                co.append(dim_of[g0 + b - 1])
                vo.append(d)

    A_ptv = sp.csr_matrix((vp, (rp, cp)), shape=(n_ptv, n))
    A_oar = sp.csr_matrix((vo, (ro, co)), shape=(n_oar, n))
    I_ptv = sp.identity(n_ptv, format="csr")
    I_oar = sp.identity(n_oar, format="csr")

    # Rows: [PTV floor | OAR ceiling | PTV ceiling]
    A = sp.bmat([
        [A_ptv, I_ptv, None, None],          # A_ptv x + u >= dmin
        [A_oar, None, -I_oar, None],         # A_oar x - v <= dmax
        [A_ptv, None, None, -I_ptv],         # A_ptv x - w <= dmax_ptv
    ], format="csc")
    nvar = n + n_ptv + n_oar + n_ptv
    lo = np.concatenate([np.full(n_ptv, DMIN), np.full(n_oar, -np.inf), np.full(n_ptv, -np.inf)])
    hi = np.concatenate([np.full(n_ptv, np.inf), np.full(n_oar, DMAX_OAR), np.full(n_ptv, DMAX_PTV)])

    # Bounds: x in [0, M], u, v, w >= 0
    lb = np.zeros(nvar)
    ub = np.full(nvar, np.inf)
    ub[:n] = MAX_INTENSITY
    A_bounds = sp.identity(nvar, format="csc")

    organ_of_oar = np.array([0] * inst.n_rows[OAR_ORGANS[0]] + [1] * inst.n_rows[OAR_ORGANS[1]])
    # normalize=True: each term is the mean squared slack of its structure (1/n).
    sz = np.array([inst.n_rows[OAR_ORGANS[0]], inst.n_rows[OAR_ORGANS[1]]], dtype=float)
    s_ptv = 1.0 / n_ptv if normalize else 1.0
    s_org = (1.0 / sz)[organ_of_oar] if normalize else np.ones(n_oar)
    w_org_vox = np.asarray(w_over_organ, dtype=float)[organ_of_oar] * s_org
    w_diag = np.concatenate([
        np.zeros(n),
        np.full(n_ptv, 2.0 * w_under * s_ptv),
        2.0 * w_org_vox,
        np.full(n_ptv, 2.0 * w_ptv_over * s_ptv),
    ])
    P = sp.diags(w_diag, format="csc")
    q = np.zeros(nvar)

    AA = sp.vstack([A, A_bounds], format="csc")
    l = np.concatenate([lo, lb])
    u_ = np.concatenate([hi, ub])

    prob = osqp.OSQP()
    prob.setup(P=P, q=q, A=AA, l=l, u=u_, eps_abs=eps, eps_rel=eps,
               max_iter=max_iter, polish=True, verbose=False)
    res = prob.solve()
    x = res.x
    if x is None or res.info.status not in ("solved", "solved inaccurate"):
        return {"status": res.info.status, "objective": float("inf")}

    xs = x[:n]
    u = x[n:n + n_ptv]
    v = x[n + n_ptv:n + n_ptv + n_oar]
    w = x[n + n_ptv + n_oar:]
    obj = float(w_under * s_ptv * np.sum(u ** 2) + np.sum(w_org_vox * v ** 2)
                + w_ptv_over * s_ptv * np.sum(w ** 2))
    # Unweighted per-structure mean squared slack: the objectives of the Pareto sweep.
    f_ptv = float(np.mean(u ** 2))
    f_bladder = float(np.mean(v[organ_of_oar == 0] ** 2))
    f_rectum = float(np.mean(v[organ_of_oar == 1] ** 2))
    D_ptv = A_ptv @ xs
    D_oar = A_oar @ xs
    return {"status": res.info.status, "objective": obj, "f_ptv": f_ptv, "f_bladder": f_bladder,
            "f_rectum": f_rectum, "x": xs, "u": u, "v": v, "w": w,
            "D_ptv": D_ptv, "D_oar": D_oar, "n_ptv": n_ptv, "n_oar": n_oar,
            "organ_of_oar": organ_of_oar, "dimlets": sorted(dim_of, key=dim_of.get)}
