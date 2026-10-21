#ifndef IMRT_FMO_H
#define IMRT_FMO_H

#include "imrt_fmo_source.h"
#include <vector>
#include <utility>
#include <memory>
#include <string>

namespace ampl { class AMPL; }

namespace emili {
namespace imrt {

/**
 * ImrtFmoSolver
 *
 * Solves the FMO QP for a set of active angles via AMPL + Gurobi, reading
 * the shared model at ampl_gurobi/fmo.mod (kept as the single source of
 * truth for the objective/constraints instead of duplicating it here):
 *
 *   min   w_under * ||u||^2  +  w_over * ||v||^2  [+  w_ptv_over * ||w||^2]
 *   s.t.  D_ptv * x + u  >=  Dmin          (PTV underdose slack)
 *         D_oar * x - v  <=  Dmax          (OAR overdose slack)
 *         D_ptv * x - w  <=  1.07*Dmin     (PTV overdose slack, if w_ptv_over>0)
 *         u, v, w  >=  0
 *         0  <=  x_j  <=  M
 *
 * The QP is built at each solve() call using only the K active angles'
 * beamlets (K x per_angle variables instead of n_candidates x per_angle).
 *
 * Optional per-organ model: setting EMILI_FMO_MODEL=organs loads
 * ampl_gurobi/fmo_organs.mod instead (one overdose weight per OAR organ, from
 * EMILI_W_OVER_ORGANS, comma-separated in oar_bounds_ order; defaults to
 * w_over for every organ; EMILI_FMO_NORMALIZE=1 divides each term by its
 * voxel count). Unset or "baseline" keeps fmo.mod and the default behavior.
 *
 * Depends on IFmoDataSource rather than a concrete instance format, so this
 * solver stays agnostic to on-disk layout details such as CerrFmoSource's
 * variable beamlet count per angle.
 *
 * Keeps a persistent ampl::AMPL member (environment + model read once in
 * the constructor) since BAO calls solve() potentially thousands of times;
 * only the per-solve data (active dimlets, sparse dose sets) is refreshed.
 */
/**
 * Resultado de un solve(): intensidades + objetivo agregado + el desglose
 * por órgano que el objetivo agregado esconde (u_b/v_b, ver fmo.mod). Cada
 * entrada de ptv_underdose_sq/oar_overdose_sq es sum(u_b^2) / sum(v_b^2)
 * sobre los boxets de ESE órgano — mismas unidades y mismo orden que
 * ptvOrganNames()/oarOrganNames(), así que zip(nombre, valor) da el
 * desglose completo sin ambigüedad.
 */
struct FmoResult {
    std::vector<double> intensities;
    double objective;
    std::vector<double> ptv_underdose_sq;  // por órgano PTV, mismo orden que ptvOrganNames()
    std::vector<double> oar_overdose_sq;   // por órgano OAR, mismo orden que oarOrganNames()
};

class ImrtFmoSolver {
public:
    explicit ImrtFmoSolver(const IFmoDataSource& source);
    ~ImrtFmoSolver();

    // Solve FMO for the given active angle indices (0-based, into the BAO
    // angle catalog -- same indexing convention the data source expects).
    FmoResult solve(const std::vector<int>& active_angles);

    int nBeamlets() const { return source_.n_dimlets(); }
    bool isReady()  const { return ready_; }

    // Nombres de órganos en el mismo orden que los vectores de FmoResult —
    // constantes durante toda la corrida (fijados en precompute()).
    std::vector<std::string> ptvOrganNames() const;
    std::vector<std::string> oarOrganNames() const;

    // When on, solve() prints per-organ boxet/dose-entry counts before each
    // real Gurobi call (skipped entirely on cache hits upstream in
    // BaoProblem, since those never reach solve()) -- a sanity check the
    // project's advisor asked for explicitly: confirms each solve is
    // actually built from fresh per-configuration data, not silently
    // reusing the same matrix across different active-angle sets.
    void setVerbose(bool v) { verbose_ = v; }

private:
    struct OrganBounds { std::string name; int row_off; int n_boxets; };

    const IFmoDataSource& source_;
    bool ready_;
    bool verbose_ = false;

    int n_ptv_, n_oar_;
    std::vector<double> dmin_;      // per PTV boxet row, length n_ptv_
    std::vector<double> dmax_;      // per OAR boxet row, length n_oar_
    std::vector<double> dmax_ptv_;  // per PTV boxet row, length n_ptv_
    std::vector<OrganBounds> ptv_bounds_;  // per-organ row ranges within dmin_/dmax_ptv_
    std::vector<OrganBounds> oar_bounds_;  // per-organ row ranges within dmax_

    std::unique_ptr<ampl::AMPL> ampl_;

    void precompute();
    void initAmpl();
    void printSolveDims(int n_active_beamlets,
                         const std::vector<int>& ptv_nnz,
                         const std::vector<int>& oar_nnz) const;
};

} // namespace imrt
} // namespace emili

#endif // IMRT_FMO_H
