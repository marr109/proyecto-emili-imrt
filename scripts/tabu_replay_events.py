#!/usr/bin/env python3
"""Rebuild the tabu-memory events of finished `tabubao` runs from their trajectory.csv.

trajectory.csv lists every evaluated candidate in order (angles, objective, cached), and the
angle-shift neighbourhood is deterministic (slots 0..K-1, -step then +step, collisions skipped),
so the search can be replayed offline using the logged objectives as an oracle. The replay is
validated row by row against the logged evaluations; a run that does not match is reported and
written with the suffix `.MISMATCH`.

Writes `tabu_events.csv` next to each trajectory.csv:
  tabu_attr runs (attribute memory):   banned angles before/after, inserted, expired
  tabu runs (configuration memory):    stored configurations, tenure, inserted, expired

Usage: tabu_replay_events.py [experiments/new]
"""
import csv
import os
import re
import sys


def read_traj(path):
    rows = []
    for r in csv.DictReader(open(path)):
        rows.append((tuple(int(float(a)) for a in r["angles_deg"].split(";")), float(r["objective"]), r["cached"] == "true"))
    return rows


class Oracle:
    """Replays the logged evaluations in order and checks them against the simulated search."""

    def __init__(self, rows):
        self.rows, self.ptr, self.real, self.table = rows, 0, 0, {}
        for a, o, _ in rows:
            self.table[tuple(sorted(a))] = o
        self.seen = set()

    def evaluate(self, angles):
        """Returns objective, or None when the log ends. Raises on a mismatch."""
        if self.ptr >= len(self.rows):
            return None
        a, o, cached = self.rows[self.ptr]
        if tuple(angles) != a:
            raise ValueError("eval %d: simulated %s, logged %s" % (self.ptr + 1, angles, a))
        key = tuple(sorted(angles))
        if (key in self.seen) != cached:
            raise ValueError("eval %d: cached flag differs" % (self.ptr + 1))
        self.seen.add(key)
        self.ptr += 1
        if not cached:
            self.real += 1
        return o


def neighbours(cur, step):
    for slot in range(len(cur)):
        for d in (-step, step):
            c = (cur[slot] + d) % 360
            if c in cur:
                continue
            v = list(cur)
            v[slot] = c
            yield tuple(v)


def fmt(v):
    return ";".join(str(x) for x in v)


class AttrMemory:
    def __init__(self, tenure):
        self.tenure, self.removed, self.cur = tenure, [], None

    def banned(self):
        return list(self.removed)

    def is_tabu(self, cand):
        return any(a in self.removed and a not in self.cur for a in cand)

    def forbid(self, s):
        if self.cur is not None:
            for a in self.cur:
                if a not in s:
                    self.removed.append(a)
                    if len(self.removed) > self.tenure:
                        self.removed.pop(0)
        self.cur = tuple(s)


class ConfigMemory:
    """Mirror of AdaptiveBaoTabuMemory: circular array of slot-ordered configurations.

    The array is replicated as is, including how it is re-laid-out when the tenure changes, because
    the entries that isTabu() scans (and so what counts as banned) depend on it.
    """

    def __init__(self, tmin=3, tmax=8):
        t0 = (tmin + tmax) // 2   # AdaptiveBaoTabuMemory starts at the midpoint, not at tmin
        self.tenure, self.tmin, self.tmax = t0, tmin, tmax
        self.mem, self.head, self.count, self.quiet = [None] * t0, 0, 0, 0

    def _scan(self):
        return [self.mem[(self.head - 1 - i + self.tenure) % self.tenure] for i in range(self.count)]

    def banned(self):
        return [fmt(c) for c in self._scan() if c is not None]

    def is_tabu(self, cand):
        return tuple(cand) in [c for c in self._scan() if c is not None]

    def forbid(self, s):
        if self.is_tabu(s):
            if self.tenure < self.tmax:
                self._resize(self.tenure + 1)
            self.quiet = 0
        else:
            self.quiet += 1
            if self.quiet >= 2 * self.tenure and self.tenure > self.tmin:
                self._resize(self.tenure - 1)
                self.quiet = 0
        self.mem[self.head] = tuple(s)
        self.head = (self.head + 1) % self.tenure
        if self.count < self.tenure:
            self.count += 1

    def _resize(self, n):
        new = [None] * n
        entries = min(self.count, n)
        for i in range(entries):
            new[(n - 1 - i + n) % n] = self.mem[(self.head - 1 - i + self.tenure) % self.tenure]
        self.mem, self.head, self.count, self.tenure = new, entries % n, entries, n


def replay(rows, step, first, mem, attr):
    orc = Oracle(rows)
    cur = rows[0][0]
    cur_obj = orc.evaluate(cur)
    best = cur_obj
    ev = []

    def snap():
        return mem.banned()

    def log(it, kind, angles, obj, before, after):
        before_s, after_s = set(map(str, before)), set(map(str, after))
        ins = [x for x in after if str(x) not in before_s]
        exp = [x for x in before if str(x) not in after_s]
        ev.append([it, kind, fmt(angles), round(obj, 6), "|".join(map(str, before)) if not attr else fmt(before),
                   "|".join(map(str, after)) if not attr else fmt(after),
                   "|".join(map(str, ins)) if not attr else fmt(ins), "|".join(map(str, exp)) if not attr else fmt(exp),
                   mem.tenure, orc.real])

    before = snap()
    mem.forbid(cur)
    log(0, "start", cur, cur_obj, before, snap())
    it = 0
    while True:
        chosen = best_any = None
        stop = False
        for cand in neighbours(cur, step):
            o = orc.evaluate(cand)
            if o is None:
                stop = True
                break
            asp = o < best
            if asp:
                best = o
            if best_any is None or o < best_any[1]:
                best_any = (cand, o)
            if not (asp or not mem.is_tabu(cand)):
                b = snap()
                log(it + 1, "tabu_skip", cand, o, b, b)
                continue
            if chosen is None or o < chosen[1]:
                chosen = (cand, o)
            if first and o < cur_obj:
                break
        if stop:
            break
        if chosen is None:
            chosen = best_any
        if chosen is None:
            break
        cur, cur_obj = chosen
        before = snap()
        mem.forbid(cur)
        it += 1
        log(it, "move", cur, cur_obj, before, snap())
    return ev, orc


HEADER = ["iter", "event", "angles", "objective", "banned_before", "banned_after", "inserted", "expired", "tenure", "real_solves"]


def main(root):
    jobs = []
    for d, _, files in os.walk(root):
        if "trajectory.csv" not in files:
            continue
        m = re.search(r"(tabu_attr)/nangshift(\d+)-t(\d+)/", d + "/") or re.search(r"(tabu)/nangshift(\d+)/", d + "/")
        if not m or "tabu_attr_trace" in d or "catalog_restricted" in d:
            continue
        strat = "first" if "/first/" in d + "/" else "best"
        jobs.append((d, m.group(1), int(m.group(2)), int(m.group(3)) if m.group(1) == "tabu_attr" else None, strat))
    ok = bad = 0
    for d, fam, step, tenure, strat in sorted(jobs):
        rows = read_traj(os.path.join(d, "trajectory.csv"))
        mem = AttrMemory(tenure) if fam == "tabu_attr" else ConfigMemory()
        try:
            ev, orc = replay(rows, step, strat == "first", mem, fam == "tabu_attr")
            if orc.ptr != len(rows):
                raise ValueError("replay consumed %d of %d logged evaluations" % (orc.ptr, len(rows)))
            out, ok = os.path.join(d, "tabu_events.csv"), ok + 1
        except ValueError as e:
            print("MISMATCH", d, e)
            out, bad = os.path.join(d, "tabu_events.csv.MISMATCH"), bad + 1
            ev = []
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(HEADER)
            w.writerows(ev)
    print("replayed ok: %d, mismatches: %d" % (ok, bad))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "experiments/new")
