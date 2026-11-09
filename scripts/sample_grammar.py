#!/usr/bin/env python3
"""Sample random commands from grammars/imrt_bao.xml (smoke test of the grammar).

  python scripts/sample_grammar.py [n=12] [seed=0]   -> prints one algorithm command per line
Parameters are drawn uniformly in their range; the sa constraint t_start > t_end is enforced by resampling.
"""
import random
import sys
import xml.etree.ElementTree as ET

PATH = "grammars/imrt_bao.xml"
NS = "{bao}"


def load():
    root = ET.parse(PATH).getroot()
    der = root.find("{grammar}derivations")
    return {c.tag.replace(NS, ""): c for c in der}


def draw(el, rng):
    t = el.get("type")
    lo, hi, st = float(el.get("min")), float(el.get("max")), float(el.get("stepIfEnumerated"))
    k = rng.randint(0, int(round((hi - lo) / st)))
    v = lo + k * st
    return str(int(v)) if t == "int" else ("%g" % round(v, 6))


def expand(name, rules, rng):
    el = rules[name]
    if el.get("type"):
        return draw(el, rng)
    alts, cur = [], []
    for ch in el.iter():
        pass
    parts, alts = [], []
    cur = []
    if el.text and el.text.strip():
        cur.append(el.text)
    for ch in el:
        if ch.tag == "or":
            alts.append(cur); cur = []
        else:
            cur.append(expand(ch.tag.replace(NS, ""), rules, rng))
        if ch.tail:
            cur.append(ch.tail)
    alts.append(cur)
    return "".join(rng.choice(alts))


def sample(rules, rng):
    while True:
        s = expand("main_algorithm", rules, rng)
        s = " ".join(s.split())
        toks = s.split()
        for key in ("sa_metropolis", "saacc"):
            if key in toks:
                i = toks.index(key)
                if float(toks[i + 1]) <= float(toks[i + 2]):
                    break
        else:
            return s


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    rules = load()
    for _ in range(n):
        print(sample(rules, rng))
