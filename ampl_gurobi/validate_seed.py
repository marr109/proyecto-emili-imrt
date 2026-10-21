"""Validate the CERR pipeline against a logged seed.
Evaluates every distinct angle set in the seed trajectory with EQUAL OAR weights
(0.5 / 0.5, the seed's w_over) and compares with the logged objective."""
import csv, os, time, sys
from cerr_fmo_ampl import CerrInstance
from solvers import build_and_solve

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Usage: python validate_seed.py [path/to/trajectory.csv]  (default: a local-search seed under experiments/)
SEED = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    REPO, "experiments", "local_search", "nangshift5", "restricted", "data", "first", "seed01", "trajectory.csv")

inst = CerrInstance()
rows = list(csv.DictReader(open(SEED)))
seen = {}
for r in rows:
    key = tuple(sorted(int(a) for a in r["angles_deg"].split(";")))
    if key not in seen:
        seen[key] = float(r["objective"])
print(f"distinct seed sets: {len(seen)}")
worst = 0.0
for key, logged in seen.items():
    t = time.time()
    out = build_and_solve(inst, key, [0.5, 0.5])
    rel = abs(out["objective"] - logged) / abs(logged)
    worst = max(worst, rel)
    print(f"{str(list(key)):<22} logged={logged:14.6f}  ours={out['objective']:14.6f}  rel={rel:.2e}  {out['status']}  ({time.time()-t:.1f}s)", flush=True)
print(f"worst relative error = {worst:.2e}")
