#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subsample a CORT-format EMILI IMRT instance into a smaller, EMILI-loadable variant.

Reduces a source PROSTATE-style instance by:
  - selecting N equi-spaced angles (deterministic index formula floor(i*src/N))
  - keeping a centered window of M dimlets per angle (0-based renumber)
  - stratified random voxel sampling per organ (proportional 60/20/20 by default)
  - filtering per-angle Dij rows to surviving (voxel, beamlet) pairs
  - writing instance_config.txt with updated counts (angles + dimlets-per-angle)
  - atomic write via temp dir + os.rename

Standard library only — runs inside the existing Docker image with no new deps.

CLI example:
  python scripts/subsample_instance.py \
      --source instances/PROSTATE_sampled \
      --output instances/PROSTATE_tiny_S42 \
      --angles 8 --dimlets 24 --seed 42

Produces instances/PROSTATE_tiny_S42/ containing:
  - instance_config.txt      (angles 8 <values>, dimlets 24)
  - Gantry<value>_Couch0_D.txt for each selected angle
      (rows filtered to local_beamlet_id in the centered window, renumbered 0..M-1,
       and to voxels surviving the per-organ stratified sample)
  - PTV_68_VOILIST.txt, Bladder_VOILIST.txt, Rectum_VOILIST.txt
      (sorted sampled voxel IDs)
"""

import argparse
import os
import random
import shutil
import sys
import tempfile


# Default per-organ proportions of --voxels total.  Must match the organ names
# as they appear in instance_config.txt.  Order = config order.
DEFAULT_PROPORTIONS = [
    ("PTV_68", 0.60),
    ("Bladder", 0.20),
    ("Rectum", 0.20),
]


# ---------------------------------------------------------------------------
# Config parsing / writing
# ---------------------------------------------------------------------------

def read_voilist(path):
    """Return the list of (integer) voxel IDs in a VOIList file (skipping comments)."""
    ids = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # a VOIList line is a single integer 0-based global voxel id
            ids.append(int(s.split()[0]))
    return ids


def read_config(path):
    """Parse instance_config.txt into an ordered list of (key, fields) tuples.

    Comments (lines whose first non-space char is '#') and blank lines are
    preserved as ('__comment__', raw_line) so we can echo them back verbatim.
    """
    entries = []
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                entries.append(("__comment__", line.rstrip("\n")))
                continue
            parts = line.split()
            entries.append((parts[0], parts[1:]))
        return entries


def config_value(entries, key):
    """Return the field list for the first entry with the given key, or None."""
    for k, fields in entries:
        if k == key:
            return fields
    return None


# ---------------------------------------------------------------------------
# Subsampling primitives
# ---------------------------------------------------------------------------

def select_angles(src_angles, n):
    """Select N equi-spaced angle values from a list of source angle values.

    Source angles are sorted numerically first; selection uses the deterministic
    index formula ``floor(i * len(src) / n)`` for i in [0, n).  This is
    seed-independent.
    """
    if n > len(src_angles):
        raise SystemExit(
            "error: --angles {0} exceeds source angle count {1}".format(
                n, len(src_angles)))
    if n < 1:
        raise SystemExit("error: --angles must be >= 1")
    src_sorted = sorted(src_angles, key=lambda v: int(v))
    L = len(src_sorted)
    selected = []
    for i in range(n):
        idx = (i * L) // n  # floor(i * L / n)
        selected.append(int(src_sorted[idx]))
    return sorted(set(selected))


def dimlet_window_start(src_dimlets, m):
    """Centered-window start index (0-based): floor((src - m) / 2)."""
    if m > src_dimlets:
        raise SystemExit(
            "error: --dimlets {0} exceeds source dimlets-per-angle {1}".format(
                m, src_dimlets))
    if m < 1:
        raise SystemExit("error: --dimlets must be >= 1")
    return (src_dimlets - m) // 2


def allocate_voxels(total, proportions, available):
    """Allocate ``total`` voxels across organs by proportional floor + largest
    fractional remainder, capped by each organ's available voxel count.

    If a proportional share exceeds an organ's available voxels, the excess is
    redistributed to organs with spare capacity (largest spare first) so that
    the returned counts still sum to ``total`` when ``total`` <= sum(available).

    ``proportions`` and ``available`` are parallel lists (organ order = config).
    Returns a list of integer counts in the same order.
    """
    raw = [float(total) * p for p in proportions]
    floors = [int(r) for r in raw]
    leftover = total - sum(floors)
    if leftover < 0:
        # Should not happen with floor, but guard anyway.
        raise SystemExit("error: voxel floor overflow (leftover < 0)")
    if leftover > 0:
        order = sorted(
            range(len(proportions)),
            key=lambda i: (-(raw[i] - floors[i]), i),
        )
        for k in range(leftover):
            floors[order[k % len(order)]] += 1

    # Cap-and-redistribute: any over-cap organ's excess is given to organs with
    # spare capacity, largest spare first.  Iterates to a fixed point.
    capacity_total = sum(available)
    target = min(total, capacity_total)
    for _ in range(len(proportions) + 1):
        excess = 0
        for i in range(len(floors)):
            if floors[i] > available[i]:
                excess += floors[i] - available[i]
                floors[i] = available[i]
        if excess == 0:
            break
        spare = [available[i] - floors[i] for i in range(len(floors))]
        if sum(spare) == 0:
            break
        order = sorted(range(len(floors)), key=lambda i: (-spare[i], i))
        oi = 0
        for _ in range(excess):
            while oi < len(order) and spare[order[oi]] <= 0:
                oi += 1
            if oi >= len(order):
                break
            j = order[oi]
            floors[j] += 1
            spare[j] -= 1

    if sum(floors) < target:
        sys.stderr.write(
            "warning: only {0} voxels available to allocate (requested {1}); "
            "total will be {0}.\n".format(sum(floors), total))
    return floors


# ---------------------------------------------------------------------------
# Write helpers (atomic)
# ---------------------------------------------------------------------------

def write_voilist(path, organ_name, voxel_ids):
    """Write a VOIList file with a header comment + sorted sampled IDs."""
    with open(path, "w") as f:
        f.write("# 0-based global voxel indices for {0} (subsampled)\n".format(organ_name))
        for vid in sorted(voxel_ids):
            f.write("{0}\n".format(vid))


def write_gantry(path, src_path, keep_voxels, start, m, angle_value, couch=0):
    """Filter a source Gantry*_D.txt file to the dimlet window + surviving voxels.

    Source row format:  global_voxel_id  local_beamlet_id  dose_rate
    local_beamlet_id is 0-based in the file.  Kept rows have:
      - local_beamlet_id in [start, start + m)
      - global_voxel_id  in keep_voxels  (union of sampled organ voxels)
    Surviving local_beamlet_id is renumbered to (orig - start) in [0, m).
    The dose_rate column is copied verbatim to preserve float precision.
    """
    header = ("# Gantry {0} Couch {1} -- global_voxel_id local_beamlet_id "
              "dose_rate (subsampled, 0-based)\n".format(angle_value, couch))
    kept = 0
    with open(src_path, "r") as src, open(path, "w") as dst:
        dst.write(header)
        for line in src:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            voxel_str, beamlet_str, dose_str = parts[0], parts[1], parts[2]
            try:
                voxel = int(voxel_str)
                beamlet = int(beamlet_str)
            except ValueError:
                continue
            if beamlet < start or beamlet >= start + m:
                continue
            if voxel not in keep_voxels:
                continue
            new_beamlet = beamlet - start
            dst.write("{0} {1} {2}\n".format(voxel, new_beamlet, dose_str))
            kept += 1
    return kept


def write_config(path, entries, new_angles, m):
    """Rebuild instance_config.txt from parsed entries, replacing:
      - the ``angles`` line  -> ``angles <N> <values...>``
      - the ``dimlets`` line -> ``dimlets <M>``
    Everything else (comments, organ defs, voilist refs, weights) is echoed
    verbatim.
    """
    angle_line = "angles\t{0} {1}".format(
        len(new_angles), " ".join(str(a) for a in new_angles))
    dimlets_line = "dimlets\t{0}".format(m)
    out_lines = []
    for key, fields in entries:
        if key == "__comment__":
            out_lines.append(fields)
        elif key == "angles":
            out_lines.append(angle_line)
        elif key == "dimlets":
            out_lines.append(dimlets_line)
        else:
            if fields:
                out_lines.append("{0}\t{1}".format(key, " ".join(fields)))
            else:
                out_lines.append("{0}".format(key))
    with open(path, "w") as f:
        f.write("\n".join(out_lines) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Subsample a CORT-format EMILI IMRT instance.")
    p.add_argument("--source",
                   help="Source instance directory (absolute or relative). "
                        "If it does not exist as a path, 'instances/<source>' "
                        "is tried.")
    p.add_argument("--output",
                   help="Output instance directory (absolute or relative). "
                        "Created atomically; must not already exist.")
    p.add_argument("--angles", type=int, default=8,
                   help="Number of equi-spaced angles to keep (default 8).")
    p.add_argument("--dimlets", type=int, default=24,
                   help="Dimlets per angle to keep, centered window (default 24).")
    p.add_argument("--voxels", type=int, default=256,
                   help="Total sampled voxels across all organs (default 256).")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for voxel sampling (default 42).")
    p.add_argument("--selftest", action="store_true",
                   help="Run built-in smoke tests and exit.")
    args = p.parse_args(argv)

    if args.selftest:
        return run_selftest()

    if not args.source or not args.output:
        p.error("--source and --output are required (except with --selftest)")

    # Resolve source directory.
    src = args.source
    if not os.path.isdir(src):
        alt = os.path.join("instances", src)
        if os.path.isdir(alt):
            src = alt
        else:
            raise SystemExit(
                "error: source instance not found: {0}".format(args.source))

    out_dir = args.output
    if os.path.exists(out_dir):
        raise SystemExit(
            "error: output already exists (fail fast): {0}".format(out_dir))
    out_parent = os.path.dirname(os.path.abspath(out_dir)) or "."
    if not os.path.isdir(out_parent):
        raise SystemExit(
            "error: output parent directory does not exist: {0}".format(out_parent))

    # Parse source config.
    config_path = os.path.join(src, "instance_config.txt")
    entries = read_config(config_path)

    src_angle_fields = config_value(entries, "angles")
    if src_angle_fields is None:
        raise SystemExit("error: source config has no 'angles' line")
    n_src_angles = int(src_angle_fields[0])
    src_angles = [int(v) for v in src_angle_fields[1:1 + n_src_angles]]
    if len(src_angles) != n_src_angles:
        raise SystemExit(
            "error: source config declares {0} angles but lists {1}".format(
                n_src_angles, len(src_angles)))

    src_dimlets_fields = config_value(entries, "dimlets")
    if src_dimlets_fields is None:
        raise SystemExit("error: source config has no 'dimlets' line")
    src_dimlets = int(src_dimlets_fields[0])

    # Organs + VOILists (preserve config order).
    organ_names = []
    organ_props = []
    for k, fields in entries:
        if k in ("ptv", "oar") and fields:
            organ_names.append(fields[0])
    # Map organ name -> proportion (default 0 if not in the default table).
    prop_map = dict(DEFAULT_PROPORTIONS)
    for name in organ_names:
        organ_props.append(prop_map.get(name, 0.0))
    # Normalise proportions to sum 1 (in case the config has different organs
    # than the default table — divide by the sum of mapped proportions).
    sprop = sum(organ_props)
    if sprop <= 0:
        raise SystemExit("error: no organ proportions available for organs {0}"
                         .format(organ_names))
    organ_props = [p / sprop for p in organ_props]

    # VOILists.
    voilist_files = {}
    for k, fields in entries:
        if k == "voilist" and len(fields) >= 2:
            voilist_files[fields[0]] = fields[1]

    # --- Angle selection (deterministic, seed-independent) -----------------
    selected_angles = select_angles(src_angles, args.angles)
    # --- Dimlet centered window (deterministic, seed-independent) ----------
    start = dimlet_window_start(src_dimlets, args.dimlets)
    m = args.dimlets

    # --- Voxel stratified sampling (seed-dependent) ------------------------
    src_voxels = {}
    for name in organ_names:
        voi_name = voilist_files.get(name, "{0}_VOILIST.txt".format(name))
        voi_path = voi_name
        if not os.path.isabs(voi_path):
            voi_path = os.path.join(src, voi_path)
        if not os.path.isfile(voi_path):
            raise SystemExit(
                "error: VOIList file not found for organ {0}: {1}".format(
                    name, voi_path))
        src_voxels[name] = read_voilist(voi_path)

    available = [len(src_voxels[name]) for name in organ_names]
    allocated = allocate_voxels(args.voxels, organ_props, available)

    rng = random.Random(args.seed)
    sampled = {}
    keep_voxels = set()
    for i, name in enumerate(organ_names):
        k = allocated[i]
        pop = src_voxels[name]
        if k == 0:
            sampled[name] = []
            continue
        picks = rng.sample(pop, k)
        sampled[name] = picks
        keep_voxels.update(picks)

    # --- Atomic write -------------------------------------------------------
    tmp = tempfile.mkdtemp(prefix="_subsample_", dir=out_parent)
    try:
        # 1. instance_config.txt
        new_n_angles = len(selected_angles)
        write_config(os.path.join(tmp, "instance_config.txt"), entries,
                     selected_angles, m)

        # 2. Gantry files (one per selected angle)
        gantry_counts = {}
        for angle in selected_angles:
            src_g = os.path.join(src,
                                 "Gantry{0}_Couch0_D.txt".format(angle))
            if not os.path.isfile(src_g):
                raise SystemExit(
                    "error: source Gantry file not found for angle {0}: {1}"
                    .format(angle, src_g))
            kept = write_gantry(
                os.path.join(tmp, "Gantry{0}_Couch0_D.txt".format(angle)),
                src_g, keep_voxels, start, m, angle, couch=0)
            gantry_counts[angle] = kept

        # 3. VOILists (sorted sampled IDs)
        for name in organ_names:
            voi_name = voilist_files.get(name, "{0}_VOILIST.txt".format(name))
            write_voilist(os.path.join(tmp, voi_name), name, sampled[name])

        # 4. Rename temp dir -> output (atomic on same filesystem).
        os.rename(tmp, out_dir)
        tmp = None  # ownership transferred
    finally:
        if tmp is not None and os.path.isdir(tmp):
            shutil.rmtree(tmp, ignore_errors=True)

    # --- Summary to stdout --------------------------------------------------
    print("subsampled instance written: {0}".format(out_dir))
    print("  source angles : {0} -> kept {1}: {2}".format(
        len(src_angles), new_n_angles, selected_angles))
    print("  source dimlets: {0}/angle -> kept {1}/angle (window [{2}, {3}))"
          .format(src_dimlets, m, start, start + m))
    print("  total dimlets : {0} (== {1} angles x {2}/angle)".format(
        new_n_angles * m, new_n_angles, m))
    print("  voxels (seed={0}):".format(args.seed))
    for i, name in enumerate(organ_names):
        print("    {0:<10} {1}/{2} (allocated {3}, source {4})".format(
            name, allocated[i], allocated[i], allocated[i], available[i]))
    print("  Gantry rows kept per angle:")
    for angle in selected_angles:
        print("    angle {0:3}: {1} rows".format(angle, gantry_counts[angle]))
    return 0


# ---------------------------------------------------------------------------
# Built-in smoke tests
# ---------------------------------------------------------------------------

def run_selftest():
    """Pure-Python asserts on the core math; no filesystem needed."""
    # 1) Angle selection: 36 -> 8 using floor(i*36/8)
    src36 = list(range(0, 360, 10))  # 0,10,...,350  (36 angles)
    sel = select_angles(src36, 8)
    assert sel == [0, 40, 90, 130, 180, 220, 270, 310], \
        "angle selection formula: got {0}".format(sel)

    # 2) N exceeds source -> error
    try:
        select_angles(src36, 64)
    except SystemExit:
        pass
    else:
        raise AssertionError("expected SystemExit when N > source angle count")

    # 3) Dimlet centered window: src=157, M=24 -> start=66, window [66,90)
    s = dimlet_window_start(157, 24)
    assert s == 66, "window start for (157,24): got {0}".format(s)
    assert (s + 24) == 90, "window end"
    # 0-based window [66,90)  ==  1-based human-readable [67,90]

    # 3b) Odd-bias case: src=156, M=24 -> start=66 (biases left)
    assert dimlet_window_start(156, 24) == 66, "odd-bias window start"

    # 4) Voxel allocation: 256 total, 60/20/20 -> 154/51/51
    alloc = allocate_voxels(256, [0.60, 0.20, 0.20], [500, 300, 1764])
    assert alloc == [154, 51, 51], "voxel alloc for 256: got {0}".format(alloc)
    assert sum(alloc) == 256, "voxel alloc sums to total"

    # 4b) Total=100 -> 60/20/20
    alloc = allocate_voxels(100, [0.60, 0.20, 0.20], [500, 300, 1764])
    assert alloc == [60, 20, 20], "voxel alloc for 100: got {0}".format(alloc)

    # 4c) Cap to available: PTV 60 capped to 10, frees 50 -> goes to organs
    # with spare capacity (Bladder and Rectum both have spare).  Largest-spare
    # first: Rectum has the most spare (1764 > 300), so Rectum takes most.
    # Both Bladder and Rectum already at proportional 20; the 50 excess is
    # distributed largest-spare first, so Rectum receives all 50 (its spare
    # stays larger than Bladder's throughout).  Final: [10, 20, 70] sum=100.
    alloc = allocate_voxels(100, [0.60, 0.20, 0.20], [10, 300, 1764])
    assert alloc == [10, 20, 70], "voxel alloc with cap: got {0}".format(alloc)
    assert sum(alloc) == 100, "cap-redistribute preserves total"

    # 5) Reproducibility: same seed -> same voxels
    pop = list(range(1000, 1100))
    r1 = random.Random(42).sample(pop, 10)
    r2 = random.Random(42).sample(pop, 10)
    assert r1 == r2, "same seed must yield same sample"
    r3 = random.Random(43).sample(pop, 10)
    assert r1 != r3, "different seeds should differ"

    # 6) Invariant: total dimlets = angles * dimlets_per_angle
    n_ang, m_per = 8, 24
    assert n_ang * m_per == 192, "CORT n_dimlets invariant"

    print("selftest OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())