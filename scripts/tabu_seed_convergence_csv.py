#!/usr/bin/env python3
"""Write convergence.csv next to each trajectory.csv of the tabu and tabu_attr runs.

One row per logged evaluation: eval index, real FMO solves so far (non-cached), objective, cached flag,
and the best objective found so far.

Usage: tabu_seed_convergence_csv.py [experiments/new]
"""
import csv
import os
import re
import sys


def main(root):
    n = 0
    for d, _, files in sorted(os.walk(root)):
        if "trajectory.csv" not in files or "tabu_attr_trace" in d or "catalog_restricted" in d:
            continue
        if not re.search(r"/(tabu_attr|tabu)/nangshift", d + "/"):
            continue
        rows, best, real = [], float("inf"), 0
        for r in csv.DictReader(open(os.path.join(d, "trajectory.csv"))):
            o = float(r["objective"])
            cached = r["cached"] == "true"
            if not cached:
                real += 1
            best = min(best, o)
            rows.append([r["eval"], real, o, r["cached"], best])
        with open(os.path.join(d, "convergence.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["eval", "real_solves", "objective", "cached", "best_so_far"])
            w.writerows(rows)
        n += 1
    print("convergence.csv written:", n)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "experiments/new")
