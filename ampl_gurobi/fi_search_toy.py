"""First-improvement local search over angle sets, exact FMO per neighbor (AMPL + Gurobi).
Toy BAO: 4 candidate angles (1 beamlet each), K=2 active. Organs: bladder (OAR 0), rectum (OAR 1).
Model: /Users/.../ampl_gurobi/fmo_organs.mod (per-organ w_over)."""
import itertools
from amplpy import AMPL

import os
MOD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fmo_organs.mod")

# Dose per unit intensity, rows = voxels, columns = candidate angles c0..c3
A_PTV = [[1.0, 0.2, 0.9, 0.4],   # PTV voxel 0
         [0.8, 0.9, 0.3, 0.5]]   # PTV voxel 1
A_OAR = [[0.7, 0.05, 0.2, 0.2],  # bladder voxel (organ 0)
         [0.1, 0.5, 0.1, 0.6]]   # rectum voxel  (organ 1)
ORGAN = [0, 1]
W_UNDER, W_PTV_OVER, W_OVER_ORGAN = 1.0, 0.3, [0.5, 0.8]
DMIN, DMAX, DMAX_PTV, M = 70.0, 30.0, 74.9, 100.0

def solve(active):
    """Exact FMO for the angle set `active` (list of candidate indices)."""
    # Build a fresh data file restricted to the active beamlets (relabelled 0..K-1).
    K = len(active)
    d_ptv = [f"{p} {j} {A_PTV[p][a]}" for p in range(2) for j, a in enumerate(active)]
    d_oar = [f"{o} {j} {A_OAR[o][a]}" for o in range(2) for j, a in enumerate(active)]
    ptv_set = " ".join(f"({p},{j})" for p in range(2) for j in range(K))
    oar_set = " ".join(f"({o},{j})" for o in range(2) for j in range(K))
    dat = f"""
param n_ptv := 2; param n_oar := 2; param n_organs := 2;
set DIMLETS := {' '.join(str(j) for j in range(K))};
param max_intensity := {M}; param w_under := {W_UNDER}; param w_ptv_over := {W_PTV_OVER};
param oar_organ := 0 {ORGAN[0]} 1 {ORGAN[1]};
param w_over_organ := 0 {W_OVER_ORGAN[0]} 1 {W_OVER_ORGAN[1]};
param dmin := 0 {DMIN} 1 {DMIN};
param dmax := 0 {DMAX} 1 {DMAX};
param dmax_ptv := 0 {DMAX_PTV} 1 {DMAX_PTV};
set PTV_DOSE := {ptv_set};
set OAR_DOSE := {oar_set};
param d_ptv := {' '.join(d_ptv)};
param d_oar := {' '.join(d_oar)};
"""
    ampl = AMPL()
    ampl.option["solver"] = "gurobi"
    ampl.option["show_stats"] = 0
    ampl.read(MOD)
    import tempfile, os
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
        f.write(dat); name = f.name
    ampl.read_data(name); os.unlink(name)
    ampl.solve()
    res = ampl.solve_result
    obj = ampl.get_objective("fmo_objective").value() if res == "solved" else float("inf")
    ampl.close()
    return obj, res

def neighbors(active, n_cand=4):
    """Swap one active angle for one inactive, in a fixed deterministic order."""
    for out in active:
        for inn in range(n_cand):
            if inn not in active:
                yield sorted([a for a in active if a != out] + [inn])

def first_improvement(start):
    current = sorted(start)
    f_cur, res = solve(current)
    trace = [(list(current), f_cur, res)]
    improved = True
    while improved:
        improved = False
        for cand in neighbors(current):
            f_c, res_c = solve(cand)
            trace.append((cand, f_c, res_c))
            if f_c < f_cur - 1e-9:          # FIRST improvement: take it, restart scan
                current, f_cur = cand, f_c
                improved = True
                break
    return current, f_cur, trace

if __name__ == "__main__":
    start = [0, 1]
    best, f_best, trace = first_improvement(start)
    print(f"start={start} -> local optimum={best}  f*={f_best:.6f}  (evaluations={len(trace)})")
    for cand, f, res in trace:
        print(f"  eval {cand}  f={f:.6f}  {res}")
    # Ground truth: enumerate all C(4,2) angle sets.
    all_sets = {tuple(s): solve(list(s))[0] for s in itertools.combinations(range(4), 2)}
    g = min(all_sets, key=all_sets.get)
    print("exhaustive:", {k: round(v, 6) for k, v in all_sets.items()})
    print(f"global optimum={list(g)}  f*={all_sets[g]:.6f}  FI reached global: {tuple(best)==g}")
