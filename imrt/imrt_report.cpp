#include "imrt_report.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace emili {
namespace imrt {

// ─────────────────────────────────────────────────────────────────────────────
// Helpers — ported verbatim from imrt_instance.cpp (pre-8c001d6).
// ─────────────────────────────────────────────────────────────────────────────

static double percentileSorted(const std::vector<double>& sorted, double p)
{
    if (sorted.empty()) return 0.0;
    // p in [0, 1].  Pre-condition: input sorted ascending.
    size_t idx = (size_t)std::floor(p * (sorted.size() - 1));
    return sorted[idx];
}

static double volumeAtOrAbove(const std::vector<double>& sorted, double dose_gy)
{
    // Returns fraction (0..1) of voxels with dose >= dose_gy.
    auto it = std::lower_bound(sorted.begin(), sorted.end(), dose_gy);
    return (double)(sorted.end() - it) / (double)sorted.size();
}

// One organ's sliced+sorted dose vector, paired with its FmoOrganRef. This
// replaces the old ImrtInstance::organs single mixed list: IFmoDataSource
// keeps PTV/OAR explicitly split (ptvOrgans()/oarOrgans()), so this struct
// re-flattens them into one ordered list (PTV first, then OAR) matching the
// old reportPlan's iteration order over inst.organs.
struct OrganBlock {
    const FmoOrganRef* ref;
    std::vector<double> sorted; // ascending
};

void reportPlan(const IFmoDataSource& source,
                const std::vector<double>& x,
                double fmo_objective,
                const std::vector<int>& angles_deg,
                std::ostream& os,
                const std::string& dvh_csv_path)
{
    // ── Compute per-boxet-row doses and slice back into per-organ vectors ───
    // computeDoses() fills ptv_dose/oar_dose as flat vectors in boxet-block
    // order: row r of ptv_dose falls within the organ whose
    // [cumulative_offset, cumulative_offset + n_boxets) slice it belongs to,
    // iterating ptvOrgans() in order (same for oar_dose / oarOrgans()).
    std::vector<double> ptv_dose, oar_dose;
    source.computeDoses(x, ptv_dose, oar_dose);

    std::vector<OrganBlock> blocks;
    {
        size_t offset = 0;
        for (const FmoOrganRef& org : source.ptvOrgans()) {
            OrganBlock b;
            b.ref = &org;
            b.sorted.assign(ptv_dose.begin() + offset,
                             ptv_dose.begin() + offset + org.n_boxets);
            std::sort(b.sorted.begin(), b.sorted.end());
            blocks.push_back(std::move(b));
            offset += org.n_boxets;
        }
    }
    {
        size_t offset = 0;
        for (const FmoOrganRef& org : source.oarOrgans()) {
            OrganBlock b;
            b.ref = &org;
            b.sorted.assign(oar_dose.begin() + offset,
                             oar_dose.begin() + offset + org.n_boxets);
            std::sort(b.sorted.begin(), b.sorted.end());
            blocks.push_back(std::move(b));
            offset += org.n_boxets;
        }
    }

    // ── PTV reference (first PTV organ) for CI/HI ────────────────────────────
    int ptv_idx = -1;
    double D_Rx = 0.0;
    for (size_t o = 0; o < blocks.size(); ++o) {
        if (blocks[o].ref->is_ptv) { ptv_idx = (int)o; D_Rx = blocks[o].ref->dmin; break; }
    }

    // ── Header ───────────────────────────────────────────────────────────────
    os << "================================================================\n";
    os << "  IMRT FMO RESULT   K=" << angles_deg.size() << " beams\n";
    os << "================================================================\n";
    os << "Selected angles (deg) : [ ";
    for (size_t i = 0; i < angles_deg.size(); ++i) {
        if (i) os << ", ";
        os << angles_deg[i];
    }
    os << " ]\n";
    os << "FMO objective f*      : " << std::fixed << std::setprecision(2)
       << fmo_objective << "\n";
    if (D_Rx > 0.0)
        os << "Prescription dose Rx  : " << std::setprecision(2) << D_Rx << " Gy\n";

    // ── Per-organ dose statistics ────────────────────────────────────────────
    os << "\n----------------------------------------------------------------\n";
    os << "  Per-organ dose statistics  (Gy)\n";
    os << "----------------------------------------------------------------\n";
    os << std::left << std::setw(11) << "Structure"
       << std::right << std::setw(8) << "Voxels"
       << std::setw(8) << "Dmin"
       << std::setw(8) << "Dmean"
       << std::setw(8) << "Dmax"
       << std::setw(8) << "D95"
       << std::setw(8) << "D5"
       << std::setw(8) << "D2" << "\n";
    os << "----------------------------------------------------------------\n";

    for (const OrganBlock& b : blocks) {
        const auto& s = b.sorted;
        if (s.empty()) continue;
        double dmin  = s.front();
        double dmax  = s.back();
        double dmean = std::accumulate(s.begin(), s.end(), 0.0) / s.size();
        double d95   = percentileSorted(s, 0.05);  // dose at 95% volume
        double d5    = percentileSorted(s, 0.95);
        double d2    = percentileSorted(s, 0.98);

        os << std::left << std::setw(11) << b.ref->name
           << std::right << std::setw(8) << b.ref->n_boxets
           << std::fixed << std::setprecision(2)
           << std::setw(8) << dmin
           << std::setw(8) << dmean
           << std::setw(8) << dmax;
        if (b.ref->is_ptv) {
            os << std::setw(8) << d95;
        } else {
            os << std::setw(8) << "-";
        }
        os << std::setw(8) << d5
           << std::setw(8) << d2 << "\n";
    }
    os << "----------------------------------------------------------------\n";

    // ── DVH constraints table ────────────────────────────────────────────────
    os << "\n----------------------------------------------------------------\n";
    os << "  DVH constraints  (clinical goals)\n";
    os << "----------------------------------------------------------------\n";
    os << std::left << std::setw(11) << "Structure"
       << std::setw(16) << "Metric"
       << std::setw(12) << "Goal"
       << std::setw(14) << "Achieved"
       << "Status\n";
    os << "----------------------------------------------------------------\n";

    for (const OrganBlock& b : blocks) {
        const FmoOrganRef& org = *b.ref;
        const auto& s = b.sorted;
        if (s.empty()) continue;

        if (org.is_ptv) {
            double d95 = percentileSorted(s, 0.05);
            double d2  = percentileSorted(s, 0.98);
            double goal_d95 = 0.95 * org.dmin;
            double goal_d2  = 1.07 * org.dmin;   // ICRU-83 hot-spot tolerance

            os << std::left << std::setw(11) << org.name
               << std::setw(16) << "D95 >= 0.95*Rx"
               << std::fixed << std::setprecision(2)
               << std::setw(12) << goal_d95
               << std::setw(14) << d95
               << (d95 >= goal_d95 ? "OK" : "VIOL") << "\n";

            os << std::left << std::setw(11) << org.name
               << std::setw(16) << "D2  <= 1.07*Rx"
               << std::fixed << std::setprecision(2)
               << std::setw(12) << goal_d2
               << std::setw(14) << d2
               << (d2 <= goal_d2 ? "OK" : "VIOL") << "\n";
        } else {
            double dmean = std::accumulate(s.begin(), s.end(), 0.0) / s.size();
            double v70 = 100.0 * volumeAtOrAbove(s, 70.0);
            double goal_v70 = 25.0;

            os << std::left << std::setw(11) << org.name
               << std::setw(16) << "Dmax <= limit"
               << std::fixed << std::setprecision(2)
               << std::setw(12) << org.dmax
               << std::setw(14) << s.back()
               << (s.back() <= org.dmax ? "OK" : "VIOL") << "\n";

            os << std::left << std::setw(11) << org.name
               << std::setw(16) << "V70 <= 25 %"
               << std::fixed << std::setprecision(2)
               << std::setw(12) << goal_v70
               << std::setw(14) << v70
               << (v70 <= goal_v70 ? "OK" : "VIOL") << "\n";

            os << std::left << std::setw(11) << org.name
               << std::setw(16) << "Dmean (info)"
               << std::setw(12) << "-"
               << std::setw(14) << dmean
               << "--" << "\n";
        }
    }
    os << "----------------------------------------------------------------\n";

    // ── Plan quality indices (CI, HI, coverage) ─────────────────────────────
    if (ptv_idx >= 0) {
        const auto& sp = blocks[ptv_idx].sorted;
        double d2  = percentileSorted(sp, 0.98);
        double d98 = percentileSorted(sp, 0.02);

        double v_ptv_at_Rx = volumeAtOrAbove(sp, D_Rx) * sp.size();

        // V_Rx = total voxels (any organ, PTV + OAR) receiving >= D_Rx.
        // Direct port of the old semantics ("all organs" = inst.organs, mixed
        // PTV/OAR): here that means summing over every block in `blocks`,
        // i.e. ptvOrgans() + oarOrgans() combined — the explicit interface
        // split doesn't change the meaning, just how the list is assembled.
        double v_total_at_Rx = 0.0;
        for (const OrganBlock& b : blocks)
            v_total_at_Rx += volumeAtOrAbove(b.sorted, D_Rx) * b.sorted.size();

        double CI = (v_ptv_at_Rx > 0.0) ? (v_total_at_Rx / v_ptv_at_Rx) : 0.0;
        double HI = (D_Rx > 0.0) ? ((d2 - d98) / D_Rx) : 0.0;
        double V95pct = 100.0 * volumeAtOrAbove(sp, 0.95 * D_Rx);

        os << "\n----------------------------------------------------------------\n";
        os << "  Plan quality indices\n";
        os << "----------------------------------------------------------------\n";
        os << "Conformity Index    CI (V_Rx / V_PTV)    = "
           << std::fixed << std::setprecision(3) << CI << "\n";
        os << "Homogeneity Index   HI (D2 - D98)/D_Rx   = "
           << std::setprecision(3) << HI << "\n";
        os << "PTV coverage        V95% of Rx           = "
           << std::setprecision(2) << V95pct << " %\n";
        os << "----------------------------------------------------------------\n";
    }

    // ── DVH CSV (optional) ──────────────────────────────────────────────────
    if (!dvh_csv_path.empty()) {
        std::ofstream csv(dvh_csv_path);
        if (csv.is_open()) {
            csv << "organ,dose_gy,volume_pct\n";
            const int N = 200;  // dose grid resolution
            for (const OrganBlock& b : blocks) {
                const auto& s = b.sorted;
                if (s.empty()) continue;
                double dmax = s.back();
                if (dmax <= 0.0) dmax = 1.0;
                for (int k = 0; k <= N; ++k) {
                    double dose = (dmax * k) / N;
                    double vol  = 100.0 * volumeAtOrAbove(s, dose);
                    csv << b.ref->name << ","
                        << std::fixed << std::setprecision(3) << dose << ","
                        << std::setprecision(3) << vol << "\n";
                }
            }
            os << "DVH curves written to: " << dvh_csv_path << "\n";
        } else {
            os << "[reportPlan] Could not open DVH CSV: " << dvh_csv_path << "\n";
        }
    }
}

} // namespace imrt
} // namespace emili
