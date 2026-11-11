#!/usr/bin/env python3
"""Check that irace/parameters.irace still mirrors grammars/imrt_bao.xml.

  python scripts/check_irace_params.py   -> prints OK / DIFF / KNOWN, exit 1 on any unexplained DIFF

parameters.irace is the hand-written record of the executed irace run (irace/parameters_g2c.txt is generated from the grammar
by scripts/grammar2code_irace.py), so this guards against drift between the two. Three checks:
  1. numeric bounds of every parameter,
  2. categorical values (step, strategy, scan order, neighbourhoods, perturbations, stepwise),
  3. dependencies: for every algorithm, the set of parameter combinations that can be active at the same time
     (grammar: expanded from the derivation structure; irace: enumerated from the `|` conditions) must be equal.
KNOWN lists deviations that are intentional, or that the executed experiment already contains.
"""
import itertools
import re
import sys
import xml.etree.ElementTree as ET

XML = "grammars/imrt_bao.xml"
IRACE = "irace/parameters.irace"
NS = "{bao}"

# irace name -> (XML name, factor applied to the irace bounds to get XML units)
NUMERIC = {
    "pert_angles": ("pert_angles", 1), "mag_angles": ("mag_angles", 1), "stagnation": ("stagnation", 1),
    "shake_levels": ("shake_levels", 1), "sample_size": ("sample_size", 1), "window_deg": ("window_deg", 5),
    "tenure_min": ("tenure_min", 1), "tenure_max": ("tenure_max", 1), "t_start": ("t_start", 1),
    "t_end": ("t_end", 1), "beta": ("beta", 1), "alpha": ("alpha", 1),
}
# parameters whose presence defines the dependency structure (same names on both sides)
GRAMMAR_ONLY_STEPS = {"1", "20"}
# alternatives the grammar offers but the executed irace run never explored (skipped in the comparison)
GRAMMAR_ONLY_PREFIXES = ("TBao_fixed", "TBao_attr", "pgreedy")
TRACK = {"scan_order", "stepwise", "pert_angles", "mag_angles", "stagnation", "shake_levels", "sample_size",
         "window_deg", "tenure_min", "tenure_max", "t_start", "t_end", "beta", "alpha"}
ALGORITHMS = {"ls": "local_search", "ils": "ils", "sa": "sa", "vnd": "vnd", "tabu": "tabu", "vns": "vns"}
KNOWN = {
    "shake_levels": "the executed irace run used (1, 4); the grammar was later capped at K-1 = 3 "
                    "because level 4 replaces all K = 4 angles",
}
NOTES = [
    "irace shares step, sample_size and window_deg across all neighbourhood/perturbation slots; the grammar draws them per use",
    "irace declares `step` unconditionally (it is ignored by vnd and by a VND-based vns), the grammar only where it is used",
    "the grammar also exposes `replacement` (distinct | replace) for the shift perturbations; the executed irace run always used distinct slots",
    "the grammar also offers the tabu memories TBao_fixed and TBao_attr (irace used only TBao_adaptive), the perturbations pgreedy and pgreedysample, and the unrestricted catalog (`irandomk unrestricted`)",
    "irace forbids nb1 == nb2 and t_start <= t_end; the grammar allows repeated neighbourhoods (t_start > t_end is left to the caller)",
]


# ---------------------------------------------------------------- grammar side
def load_xml():
    der = ET.parse(XML).getroot().find("{grammar}derivations")
    return {c.tag.replace(NS, ""): c for c in der}


def alternatives(el):
    """Split a derivation into alternatives: (list of child refs, plain text) per alternative."""
    alts = [([], [el.text or ""])]
    for ch in el:
        if ch.tag == "or":
            alts.append(([], []))
        else:
            ref = ch.tag.replace(NS, "")
            alts[-1][0].append(ref)
            alts[-1][1].append("<%s>" % ref)
        alts[-1][1].append(ch.tail or "")
    return [(refs, " ".join("".join(txt).split())) for refs, txt in alts]


def expand(tag, D):
    """Set of frozensets of TRACK parameters that can be active together under this derivation."""
    el = D[tag]
    if el.get("type") in ("int", "real"):
        return {frozenset([tag])} if tag in TRACK else {frozenset()}
    out = set()
    for refs, text in alternatives(el):
        if text.startswith(GRAMMAR_ONLY_PREFIXES):
            continue
        combos = {frozenset()}
        for r in refs:
            combos = {a | b for a in combos for b in expand(r, D)}
        out |= combos
    if tag in TRACK:  # scan_order / stepwise: presence is what matters, not the chosen value
        out = {frozenset([tag])}
    return out


# ------------------------------------------------------------------- irace side
def parse_irace():
    params, forbidden = {}, []
    in_forbidden = False
    for line in open(IRACE):
        line = line.rstrip("\n")
        if line.strip() == "[forbidden]":
            in_forbidden = True
            continue
        if in_forbidden and line.strip() and not line.startswith("#"):
            forbidden.append(line.strip())
            continue
        m = re.match(r'^(\w+)\s+"[^"]*"\s+([irc])\s+\(([^)]*)\)(?:\s*\|\s*(.*))?$', line)
        if m:
            vals = [v.strip() for v in m.group(3).split(",")]
            params[m.group(1)] = {"type": m.group(2), "vals": vals, "cond": (m.group(4) or "").strip()}
    return params, forbidden


def r2py(cond):
    cond = re.sub(r"\bc\(([^)]*)\)", r"(\1,)", cond)
    return cond.replace("%in%", " in ").replace("&", " and ").replace("|", " or ")


class Env(dict):
    def __missing__(self, key):
        return None


def irace_signatures(params):
    names = list(params)
    deps = {n: {d for d in names if d != n and re.search(r"\b%s\b" % d, params[n]["cond"])} for n in names}
    order, seen = [], set()

    def visit(n):
        if n in seen:
            return
        seen.add(n)
        for d in sorted(deps[n]):
            visit(d)
        order.append(n)
    for n in names:
        visit(n)
    sigs = {a: set() for a in ALGORITHMS}

    def rec(i, env):
        if i == len(order):
            sigs[env["algorithm"]].add(frozenset(n for n in TRACK if env[n] is not None))
            return
        n = order[i]
        p = params[n]
        active = (not p["cond"]) or bool(eval(r2py(p["cond"]), {"__builtins__": {}}, env))
        if not active:
            rec(i + 1, Env(env, **{n: None}))
        elif p["type"] == "c":
            for v in p["vals"]:
                rec(i + 1, Env(env, **{n: v}))
        else:
            rec(i + 1, Env(env, **{n: True}))
    rec(0, Env())
    return sigs


# ----------------------------------------------------------------------- checks
def check_bounds(D, params):
    bad = 0
    xb = {t: (float(e.get("min")), float(e.get("max"))) for t, e in D.items() if e.get("type") in ("int", "real")}
    for name, (xname, k) in NUMERIC.items():
        if name not in params or xname not in xb:
            print("DIFF  bounds %-13s missing in %s" % (name, "irace" if name not in params else "grammar")); bad += 1
            continue
        lo, hi = (float(x) * k for x in params[name]["vals"])
        if (lo, hi) == xb[xname]:
            print("OK    bounds %-13s %s" % (name, (lo, hi)))
        elif name in KNOWN:
            print("KNOWN bounds %-13s irace %s vs grammar %s -- %s" % (name, (lo, hi), xb[xname], KNOWN[name]))
        else:
            print("DIFF  bounds %-13s irace %s vs grammar %s" % (name, (lo, hi), xb[xname])); bad += 1
    return bad


def texts(D, tag):
    return [t for _, t in alternatives(D[tag])]


def check_categories(D, params):
    bad = 0
    want = {  # irace parameter -> expected set of grammar alternatives (first token / text)
        "step": ("step", {"5", "10"}, lambda t: t),
        "strategy": ("strategy", {"first", "best"}, lambda t: t),
        "scan_order": ("scan_order", {"fixed", "circular", "randorder"}, lambda t: t or "fixed"),
        "stepwise": ("stepwise", {"no", "yes"}, lambda t: "yes" if t else "no"),
    }
    for name, (tag, irace_set, norm) in want.items():
        got = {norm(t) for t in texts(D, tag)}
        extra = GRAMMAR_ONLY_STEPS if name == "step" else set()   # steps the executed irace run did not explore
        if extra and extra <= got:
            print("KNOWN categories step        grammar also admits steps %s, outside the executed irace space %s" % (sorted(extra), sorted(irace_set)))
            got -= extra
        ok = got == set(params[name]["vals"]) == irace_set
        print("%-5s categories %-11s irace %s vs grammar %s" % ("OK" if ok else "DIFF", name, sorted(params[name]["vals"]), sorted(got)))
        bad += not ok
    nb = {"nangshift 10": "shift10", "nangshift 5": "shift5", "nangswapsample": "sample",
          "nangswapwindow": "window", "nangflip 0": "flip", "nangshift 1": "shift1", "nangshift 20": "shift20"}
    got = set()
    for t in texts(D, "neighbourhood"):
        t = re.sub(r"<\w+>", "", t).strip()
        got.add(nb.get(t, "?" + t))
    if {"shift1", "shift20"} <= got:
        print("KNOWN categories nb1-nb3     grammar also admits shift1 and shift20, outside the executed irace space")
        got -= {"shift1", "shift20"}
    for n in ("nb1", "nb2", "nb3"):
        iv = set(params[n]["vals"]) - {"none"}
        ok = iv == got
        print("%-5s categories %-11s irace %s vs grammar %s" % ("OK" if ok else "DIFF", n, sorted(iv), sorted(got)))
        bad += not ok
    pert = {"prangshift": "shift", "prangshiftadaptive": "adaptive", "prangshiftmag": "mag", "prangshiftadaptivemag": "adaptivemag"}
    only = [t for t in texts(D, "perturbation") if t.startswith(GRAMMAR_ONLY_PREFIXES)]
    if only:
        print("KNOWN categories pert        grammar also admits %s, outside the executed irace space" % sorted({t.split()[0] for t in only}))
    got = {pert.get(t.split()[0], "?" + t.split()[0]) for t in texts(D, "perturbation") if t not in only}
    ok = got == set(params["pert"]["vals"])
    print("%-5s categories %-11s irace %s vs grammar %s" % ("OK" if ok else "DIFF", "pert", sorted(params["pert"]["vals"]), sorted(got)))
    bad += not ok
    return bad


def check_dependencies(D, params):
    bad = 0
    isig = irace_signatures(params)
    for a, tag in ALGORITHMS.items():
        gsig = expand(tag, D)
        if gsig == isig[a]:
            print("OK    dependencies %-4s %d active-parameter combinations agree" % (a, len(gsig)))
        else:
            bad += 1
            print("DIFF  dependencies %-4s" % a)
            for s in sorted(gsig - isig[a], key=sorted):
                print("        only in grammar:", sorted(s))
            for s in sorted(isig[a] - gsig, key=sorted):
                print("        only in irace:  ", sorted(s))
    return bad


def main():
    D = load_xml()
    params, forbidden = parse_irace()
    bad = check_bounds(D, params) + check_categories(D, params) + check_dependencies(D, params)
    print("\nKNOWN deviations (documented, not errors):")
    for n in NOTES:
        print("  -", n)
    print("  - irace [forbidden]:", "; ".join(forbidden))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
