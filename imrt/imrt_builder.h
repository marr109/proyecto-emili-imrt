#ifndef IMRT_BUILDER_H
#define IMRT_BUILDER_H

#include "../generalParser.h"
#include "imrt.h"
#include "imrt_bao.h"

namespace prs {
namespace imrt {

/**
 * ImrtBuilder — parses command-line tokens and constructs IMRT components.
 *
 * ── FMO problem (classic metaheuristic on beamlet intensities) ─────────────
 *   Token: imrt  <instance_dir>  [nactive <K>]  [verbose]
 *
 * ── BAO problem (combinatorial angle search + FMO solver) ─────────────────
 *   Token: baoimrt  <K>  <instance_dir>  [verbose]
 *   NOTE: on this branch the FMO solver (imrt_fmo.cpp) is a stub that fails
 *   at runtime — see ampl_gurobi/ for the working AMPL+Gurobi replacement.
 *
 * ── FMO initial solutions ─────────────────────────────────────────────────
 *   izero                  all intensities = 0
 *   iuniform  <x>          all intensities = x
 *   irandom   <max>        random intensities in [0, max]
 *
 * ── BAO initial solutions ─────────────────────────────────────────────────
 *   ifirstk                first K angles (0,1,...,K-1)
 *   irandomk               K randomly chosen angles
 *
 * ── FMO neighborhoods ─────────────────────────────────────────────────────
 *   nshift   <delta>       SingleBeamletShift
 *   nswap                  BeamletSwap
 *
 * ── BAO neighborhoods ─────────────────────────────────────────────────────
 *   nangswap               AngleSwapNeighborhood (swap one active ↔ inactive)
 *   nangshift <step>       AngleShiftNeighborhood (shift active angle ±step
 *                           positions along the degree-sorted catalog)
 *
 * ── Perturbations ─────────────────────────────────────────────────────────
 *   prandom  <k> <max>     RandomBeamletPerturbation
 *
 * ── Acceptance ────────────────────────────────────────────────────────────
 *   aimprove               ImrtImproveAccept
 *
 * ── Termination ───────────────────────────────────────────────────────────
 *   tmaxiter  <n>          ImrtMaxIterations / iterations on angle swaps
 *   tfeasible              ImrtFeasibleTermination
 */
class ImrtBuilder : public Builder {
public:
    ImrtBuilder(GeneralParserE& gp, TokenManager& tm)
        : Builder(gp, tm) {}

    virtual bool isCompatibleWith(char* problem_definition) override;
    virtual bool canOpenInstance(char* problem_definition) override;
    virtual emili::Problem* openInstance() override;

    virtual emili::InitialSolution*    buildInitialSolution()    override;
    virtual emili::Neighborhood*       buildNeighborhood()       override;
    virtual emili::Perturbation*       buildPerturbation()       override;
    virtual emili::Acceptance*         buildAcceptance()         override;
    virtual emili::Termination*        buildTermination()        override;
    virtual emili::TabuMemory*         buildTabuTenure()         override;
    virtual emili::Shake*              buildShake()              override;

private:
    emili::imrt::ImrtProblem* castProblem();
    emili::imrt::BaoProblem*  castBaoProblem();
    bool isBaoProblem();
};

} // namespace imrt
} // namespace prs

#endif // IMRT_BUILDER_H
