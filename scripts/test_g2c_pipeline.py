#!/usr/bin/env python3
"""End-to-end check of the grammar-driven irace pipeline (no irace needed).

  python scripts/test_g2c_pipeline.py [n=100] [seed=0] [--instance instances/CERR_Prostate] [--jobs 8]

For n random configurations of irace/parameters_g2c.txt (conditions and [forbidden] rules respected):
  1. grammar2code translates the configuration into an EMILI command line,
  2. build/emili parses and starts it (-it 2, so it only checks that the command is accepted, not its quality).
Prints how many were accepted per algorithm and the first failures. Exit code 1 if any command is rejected.
Environment: GRAMMAR2CODE = path to the grammar2code executable (default: grammar2code in the PATH).
"""
import argparse
import os
import random
import re
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G2C = os.environ.get("GRAMMAR2CODE", "grammar2code")
BAD = re.compile(r"FATAL|expecting|\[IMRT\]|rror|Segmentation|terminate called")
LINE = re.compile(r'^(\w+)\s+"([^"]*)"\s+([irc])\s+\(([^)]*)\)(?:\s*\|\s*(.*))?$')


def load(path):
    params, forbidden, in_forbidden = [], [], False
    for line in open(path):
        line = line.rstrip("\n")
        if line.strip() == "[forbidden]":
            in_forbidden = True
        elif in_forbidden and line.strip() and not line.startswith("#"):
            forbidden.append(line.strip())
        else:
            m = LINE.match(line)
            if m:
                vals = [v.strip() for v in m.group(4).split(",")]
                params.append({"name": m.group(1), "switch": m.group(2), "type": m.group(3), "vals": vals, "cond": (m.group(5) or "").strip()})
    return params, forbidden


def r2py(cond):
    cond = re.sub(r"\bc\(([^)]*)\)", r"(\1,)", cond)
    return cond.replace("%in%", " in ").replace("&", " and ").replace("|", " or ")


class Env(dict):
    def __missing__(self, key):
        return None


def sample(params, forbidden, rng):
    for _ in range(100):
        env, args = Env(), []
        for p in params:
            if p["cond"] and not eval(r2py(p["cond"]), {"__builtins__": {}}, env):
                continue
            if p["type"] == "c":
                v = int(rng.choice(p["vals"]))
            elif p["type"] == "i":
                v = rng.randint(int(p["vals"][0]), int(p["vals"][1]))
            else:
                v = round(rng.uniform(float(p["vals"][0]), float(p["vals"][1])), 4)
            env[p["name"]] = v
            args.append("%s%s" % (p["switch"], v))
        if not any(eval(f.replace("<=", "<="), {"__builtins__": {}}, env) for f in forbidden
                   if all(env[n] is not None for n in re.findall(r"\w+", f))):
            return env, args
    raise RuntimeError("could not sample a configuration that respects the forbidden rules")


def check(job):
    i, args, instance = job
    with tempfile.TemporaryDirectory(prefix="g2c_test_") as tmp:
        t = subprocess.run([G2C, os.path.join(ROOT, "grammars", "imrt_bao.xml"), "--depth=1", "--target_dir=" + tmp] + args,
                           capture_output=True, text=True)
        pfile = os.path.join(tmp, "params")
        if t.returncode != 0 or not os.path.exists(pfile):
            return i, None, "grammar2code failed: " + (t.stdout + t.stderr)[-200:]
        algo = open(pfile).read().split()
        env = dict(os.environ, EMILI_FMO_MODEL="organs", EMILI_W_OVER_ORGANS="0.5,0.8", EMILI_DMIN="68",
                   EMILI_W_PTV_OVER="0.5", EMILI_EXCLUDE_OVERLAP="1", EMILI_MAX_INTENSITY="100", EMILI_GUROBI_OPTIONS="threads=1")
        argv = [os.path.join(ROOT, "build", "emili"), os.path.join(ROOT, instance), "baoimrt", "4", "csv",
                os.path.join(tmp, "t.csv")] + algo + ["-it", "2", "rnds", "1"]
        try:
            p = subprocess.run(argv, env=env, cwd=ROOT, capture_output=True, text=True, timeout=60)
            out = p.stdout + p.stderr
        except subprocess.TimeoutExpired as e:
            out = ((e.stdout or b"") + (e.stderr or b"")).decode(errors="ignore") if isinstance(e.stdout, bytes) else (e.stdout or "") + (e.stderr or "")
        m = BAD.search(out)
        return i, " ".join(algo), (out[max(0, m.start() - 40):m.end() + 80].replace("\n", " ") if m else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", nargs="?", type=int, default=100)
    ap.add_argument("seed", nargs="?", type=int, default=0)
    ap.add_argument("--instance", default="instances/CERR_Prostate")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--params", default=os.path.join(ROOT, "irace", "parameters_g2c.txt"))
    a = ap.parse_args()
    params, forbidden = load(a.params)
    rng = random.Random(a.seed)
    jobs = []
    for i in range(a.n):
        env, args = sample(params, forbidden, rng)
        jobs.append((i, args, a.instance))
    with ThreadPoolExecutor(a.jobs) as ex:
        results = list(ex.map(check, jobs))
    ok, bad, per_alg = 0, [], Counter()
    for (i, args, _), (_, cmd, err) in zip(jobs, results):
        alg = args[0].split("=")[1]
        if err is None:
            ok += 1
            per_alg[alg] += 1
        else:
            bad.append((i, cmd, err))
    print("accepted: %d of %d   (branches of baomain_algorithm: %s)" % (ok, a.n, dict(sorted(per_alg.items()))))
    for i, cmd, err in bad[:8]:
        print("REJECTED #%d: %s\n   %s" % (i, cmd, err))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
