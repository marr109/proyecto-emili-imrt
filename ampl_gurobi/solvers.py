"""Solver dispatcher. EMILI_SOLVER=gurobi (default, AMPL+Gurobi) | osqp (no size limit, no license)."""
import os
def build_and_solve(inst, active, w, **kw):
    if os.environ.get("EMILI_SOLVER", "gurobi") == "osqp":
        from cerr_fmo_osqp import build_and_solve as f
    else:
        from cerr_fmo_gurobi import build_and_solve as f
    return f(inst, active, w, **kw)
