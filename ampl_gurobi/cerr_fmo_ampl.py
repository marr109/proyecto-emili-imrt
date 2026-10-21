"""CERR_Prostate FMO evaluation through AMPL + Gurobi with per-organ OAR weights.

Replicates imrt/cerr_instance.cpp (CerrFmoSource) so that every evaluation of
a set of K gantry angles solves the same QP the C++ pipeline would build:

  * beamletIndex.txt : angle_idx  first_global_beamlet  last_global_beamlet (1-based)
  * <ORGAN>.txt       : voxel list; the row of a voxel is its order of appearance
  * <ORGAN>_<a>.txt   : per-angle entries  voxel_id  local_beamlet  dose
  * global dimlet    = (first - 1) + (local_beamlet - 1)

Organ order and row offsets follow the C++ loader:
  PTV: PTVHD, PTVLD (dmin 65 Gy, dmax_ptv 1.07*65)
  OAR: BLADDER, RECTUM (dmax 50 Gy)
"""
import os
import tempfile
from collections import defaultdict

from amplpy import AMPL

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "fmo_organs.mod")
INSTANCE = os.path.expanduser("~/Documents/code/emili_imrt/instances/CERR_Prostate")

PTV_ORGANS = ["PTVHD", "PTVLD"]
OAR_ORGANS = ["BLADDER", "RECTUM"]
DMIN = float(os.environ.get("EMILI_DMIN", 65.0))      # PTV prescription (Gy); default = loader placeholder
DMAX_OAR = 50.0
DMAX_PTV = 1.07 * DMIN
MAX_INTENSITY = float(os.environ.get("EMILI_MAX_INTENSITY", 15000.0))   # per-beamlet cap; default = original placeholder
EXCLUDE_OVERLAP = os.environ.get("EMILI_EXCLUDE_OVERLAP", "0") == "1"   # OAR minus PTV
W_UNDER = 1.0
W_PTV_OVER = float(os.environ.get("EMILI_W_PTV_OVER", 0.0))


class CerrInstance:
    def __init__(self, path=INSTANCE):
        self.path = path
        self.beamlet_range = {}          # angle -> (first, last), 1-based inclusive
        self.voxel_row = {}              # organ -> {voxel_id: local row}
        self.n_rows = {}                 # organ -> number of voxels
        self.row_off = {}                # organ -> global row offset within its group
        self.cache = {}                  # angle -> {organ: [(voxel_id, local_beamlet, dose)]}
        self._load_index()
        self._load_voxels()

    def _load_index(self):
        with open(os.path.join(self.path, "beamletIndex.txt")) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                idx, s, e = map(int, line.split()[:3])
                self.beamlet_range[idx] = (s, e)

    def _load_voxels(self):
        ptv_ids = set()
        for group in (PTV_ORGANS, OAR_ORGANS):
            off = 0
            for organ in group:
                rows = {}
                with open(os.path.join(self.path, organ + ".txt")) as f:
                    for line in f:
                        if not line.strip() or line.startswith("#"):
                            continue
                        vid = int(line.split()[0])
                        if group is OAR_ORGANS and EXCLUDE_OVERLAP and vid in ptv_ids:
                            continue          # OAR minus PTV: overlap voxels count only as PTV
                        rows[vid] = len(rows)
                if group is PTV_ORGANS:
                    ptv_ids.update(rows)
                self.voxel_row[organ] = rows
                self.n_rows[organ] = len(rows)
                self.row_off[organ] = off
                off += len(rows)

    def angle(self, a):
        """Parse and cache all dose entries of one angle, for every organ."""
        if a in self.cache:
            return self.cache[a]
        out = {}
        for organ in PTV_ORGANS + OAR_ORGANS:
            entries = []
            fn = os.path.join(self.path, f"{organ}_{a}.txt")
            with open(fn) as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    v, b, d = line.split()[:3]
                    v = int(v)
                    if v in self.voxel_row[organ]:
                        entries.append((v, int(b), float(d)))
            out[organ] = entries
        self.cache[a] = out
        return out

    def dat_text(self, active, w_over_organ, w_under=None, w_ptv_over=None):
        w_under = W_UNDER if w_under is None else w_under
        w_ptv_over = W_PTV_OVER if w_ptv_over is None else w_ptv_over
        """AMPL data for the QP restricted to the beamlets of `active` angles."""
        dimlets = []                    # global dimlet ids of the active angles
        local_of = {}                   # (angle, local_beamlet) -> global dimlet
        for a in active:
            s, e = self.beamlet_range[a]
            for g in range(s - 1, e):
                local_of[(a, g - (s - 1) + 1)] = g
                dimlets.append(g)

        ptv_rows = self.n_rows[PTV_ORGANS[0]] + self.n_rows[PTV_ORGANS[1]]
        oar_rows = self.n_rows[OAR_ORGANS[0]] + self.n_rows[OAR_ORGANS[1]]

        d_ptv, d_oar = [], []
        ptv_pairs, oar_pairs = set(), set()
        for a in active:
            data = self.angle(a)
            g0 = self.beamlet_range[a][0] - 1
            for organ in PTV_ORGANS:
                off = self.row_off[organ]
                for v, b, d in data[organ]:
                    row = off + self.voxel_row[organ][v]
                    dim = g0 + b - 1
                    d_ptv.append(f"{row} {dim} {d}")
                    ptv_pairs.add((row, dim))
            for organ in OAR_ORGANS:
                off = self.row_off[organ]
                for v, b, d in data[organ]:
                    row = off + self.voxel_row[organ][v]
                    dim = g0 + b - 1
                    d_oar.append(f"{row} {dim} {d}")
                    oar_pairs.add((row, dim))

        oar_organ = [0] * self.n_rows[OAR_ORGANS[0]] + [1] * self.n_rows[OAR_ORGANS[1]]
        dimlet_set = " ".join(str(g) for g in dimlets)
        ptv_set = " ".join(f"({r},{j})" for r, j in sorted(ptv_pairs))
        oar_set = " ".join(f"({r},{j})" for r, j in sorted(oar_pairs))
        return "\n".join([
            f"param n_ptv := {ptv_rows}; param n_oar := {oar_rows}; param n_organs := 2;",
            f"set DIMLETS := {dimlet_set};",
            f"param max_intensity := {MAX_INTENSITY}; param w_under := {w_under};",
            f"param w_ptv_over := {w_ptv_over};",
            "param oar_organ := " + " ".join(f"{i} {o}" for i, o in enumerate(oar_organ)) + ";",
            f"param w_over_organ := 0 {w_over_organ[0]} 1 {w_over_organ[1]};",
            f"param dmin := " + " ".join(f"{i} {DMIN}" for i in range(ptv_rows)) + ";",
            f"param dmax := " + " ".join(f"{i} {DMAX_OAR}" for i in range(oar_rows)) + ";",
            f"param dmax_ptv := " + " ".join(f"{i} {DMAX_PTV}" for i in range(ptv_rows)) + ";",
            f"set PTV_DOSE := {ptv_set};",
            f"set OAR_DOSE := {oar_set};",
            "param d_ptv := " + " ".join(d_ptv) + ";",
            "param d_oar := " + " ".join(d_oar) + ";",
        ])

    def solve(self, active, w_over_organ):
        """Exact FMO for the angle set `active`. Returns a dict with objective and holguras."""
        dat = self.dat_text(sorted(active), w_over_organ)
        ampl = AMPL()
        ampl.option["solver"] = "gurobi"
        ampl.option["show_stats"] = 0
        ampl.read(MOD)
        with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
            f.write(dat)
            name = f.name
        ampl.read_data(name)
        os.unlink(name)
        ampl.solve()
        res = ampl.solve_result
        out = {"result": res, "objective": float("inf")}
        if res == "solved":
            out["objective"] = ampl.get_objective("fmo_objective").value()
            out["u"] = ampl.get_variable("u").to_dict()
            out["v"] = ampl.get_variable("v").to_dict()
            out["x"] = ampl.get_variable("x").to_dict()
        ampl.close()
        return out
