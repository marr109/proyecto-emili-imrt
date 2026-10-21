"""First-improvement local search on CERR_Prostate, K=4, displacement neighborhood.
Neighborhood (report definition, 2K moves): for each active angle in ascending order,
shift by -STEP then +STEP (mod 360); moves onto an already active angle are skipped.
Solver: OSQP on the exact QP (see cerr_fmo_osqp.py), per-organ OAR weights."""
import sys, time, json, os
TAG = os.environ.get("OUT_TAG", "")
from cerr_fmo_ampl import CerrInstance
from solvers import build_and_solve

START = [int(a) for a in os.environ.get("EMILI_START", "35,115,165,285").split(",")]   # default: logged seed01 start
STEP = 5
CONDITIONS = {"organ_weights_0.5_0.8": [0.5, 0.8], "control_equal_0.5_0.5": [0.5, 0.5]}
if os.environ.get("ONLY_ORGAN"):
    CONDITIONS = {k: v for k, v in CONDITIONS.items() if k.startswith("organ")}

def neighbors(active):
    for slot in sorted(active):
        for delta in (-STEP, +STEP):
            new = sorted([(a + delta) % 360 if a == slot else a for a in active])
            if len(set(new)) == len(new):
                yield new

def run(inst, w, label, log):
    cache = {}
    def f(angles):
        key = tuple(sorted(angles))
        if key not in cache:
            t = time.time()
            out = build_and_solve(inst, key, w)
            cache[key] = out["objective"]
            log.write(f"  [{label}] eval {len(cache):3d} {list(key)} f={out['objective']:.6f} "
                      f"{out['status']} ({time.time()-t:.0f}s)\n"); log.flush()
        return cache[key]
    current = sorted(START); f_cur = f(current)
    log.write(f"[{label}] start {current} f={f_cur:.6f}\n"); log.flush()
    improved = True
    while improved:
        improved = False
        for cand in neighbors(current):
            f_c = f(cand)
            if f_c < f_cur - 1e-9:
                log.write(f"[{label}] IMPROVE {current} -> {cand}: {f_cur:.6f} -> {f_c:.6f}\n"); log.flush()
                current, f_cur = cand, f_c; improved = True
                break
    log.write(f"[{label}] LOCAL OPTIMUM {current} f*={f_cur:.6f} evaluations={len(cache)}\n"); log.flush()
    return {"label": label, "optimum": current, "f_star": f_cur, "evaluations": len(cache)}

if __name__ == "__main__":
    inst = CerrInstance()
    results = []
    with open(f"fi_cerr_log{TAG}.txt", "w") as log:
        for label, w in CONDITIONS.items():
            results.append(run(inst, w, label, log))
    json.dump(results, open(f"fi_cerr_results{TAG}.json", "w"), indent=2)
    print("done", results)
