#ifndef IMRT_REPORT_H
#define IMRT_REPORT_H

#include <iostream>
#include <string>
#include <vector>
#include "imrt_fmo_source.h"

namespace emili {
namespace imrt {

/*---------------------------------------------------------------------------*
 * reportPlan — clinical-style summary of a fluence map.
 *
 * Ported verbatim (formulas unchanged) from the pre-CORT-removal
 * imrt_instance.cpp::reportPlan (see `git show 8c001d6^:imrt/imrt_instance.cpp`),
 * adapted to read organ doses from IFmoDataSource::computeDoses() (flat
 * ptv_dose/oar_dose vectors in boxet-block order, one block per organ in
 * ptvOrgans()/oarOrgans() order) instead of the old
 * ImrtInstance::computeOrganDoses() (native per-organ nesting).
 *
 * Given an intensity vector x (length source.n_dimlets()), prints to `os`:
 *   - selected angles, FMO objective value
 *   - per-organ Dmin, Dmean, Dmax, D95, D5, D2 (Gy) — PTV organs first
 *     (ptvOrgans() order), then OAR organs (oarOrgans() order)
 *   - DVH constraint table (PTV: D95>=0.95*Rx, D2<=1.07*Rx; OAR: Dmax<=limit,
 *     V70<=25%) flagged OK / VIOL
 *   - Conformity Index (CI), Homogeneity Index (HI), V95% coverage per
 *     ICRU-83, referenced to the first PTV organ's prescription dose (Rx)
 *
 * If dvh_csv_path is non-empty, also writes one row per organ x dose-grid
 * point ("organ,dose_gy,volume_pct", 200-point grid) usable for plotting
 * DVH curves.
 *---------------------------------------------------------------------------*/
void reportPlan(const IFmoDataSource& source,
                const std::vector<double>& x,
                double fmo_objective,
                const std::vector<int>& angles_deg,
                std::ostream& os,
                const std::string& dvh_csv_path = "");

} // namespace imrt
} // namespace emili

#endif // IMRT_REPORT_H
