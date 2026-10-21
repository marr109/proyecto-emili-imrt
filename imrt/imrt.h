#ifndef IMRT_H
#define IMRT_H

#include "../emilibase.h"
#include "imrt_fmo_source.h"
#include <memory>
#include <vector>

namespace emili {
namespace imrt {

class ImrtSolution; // forward declaration

/*---------------------------------------------------------------------------*
 *                              PROBLEM                                      *
 *---------------------------------------------------------------------------*/

class ImrtProblem : public emili::Problem {
protected:
    std::unique_ptr<IFmoDataSource> source_;
    std::vector<int>  angle_degrees_;     // catalog: angle_degrees_[i] = degree of candidate angle i
    std::vector<int>  active_angle_idxs_; // empty = all angles active
    bool              verbose_fmo_;
    double            last_best_fmo_;
    int               improve_count_;

    void printFmoStats(const ImrtSolution& sol, double val, double delta);

public:
    ImrtProblem(std::unique_ptr<IFmoDataSource> source, std::vector<int> angle_degrees)
        : source_(std::move(source)), angle_degrees_(std::move(angle_degrees)),
          verbose_fmo_(false), last_best_fmo_(1e18), improve_count_(0) {}

    /** f = w_under*sum(max(0,Dmin-dosis)^2) [PTV] + w_over*sum(max(0,dosis-Dmax)^2) [OAR]. */
    virtual double calcObjectiveFunctionValue(emili::Solution& solution) override;
    virtual double evaluateSolution(emili::Solution& solution) override;
    virtual int    problemSize() override { return source_->n_dimlets(); }

    const IFmoDataSource& getSource() const { return *source_; }
    int  nAngles() const { return (int)angle_degrees_.size(); }
    const std::vector<int>& getAngleDegrees() const { return angle_degrees_; }

    /** Use only the first k angles (0-based indices 0..k-1). */
    void setActiveAngles(int k);

    /** True if the given 0-based angle index is in the active set (or set is empty). */
    bool isAngleActive(int angle_idx) const;

    const std::vector<int>& getActiveAngles() const { return active_angle_idxs_; }

    void setVerboseFmo(bool v) { verbose_fmo_ = v; }
};


/*---------------------------------------------------------------------------*
 *                              SOLUTION                                     *
 *---------------------------------------------------------------------------*/

/**
 * ImrtSolution — a vector of beamlet intensities.
 * x[d] is the intensity of global dimlet d, d ∈ [0, n_dimlets).
 */
class ImrtSolution : public emili::Solution {
protected:
    std::vector<double> intensities_;

    virtual const void* getRawData()           const override;
    virtual void        setRawData(const void* data) override;

public:
    ImrtSolution(double value, const std::vector<double>& x)
        : emili::Solution(value), intensities_(x) {}

    explicit ImrtSolution(int n_dimlets)
        : emili::Solution(0.0), intensities_(n_dimlets, 0.0) {}

    const std::vector<double>& getIntensities() const { return intensities_; }
    void  setIntensities(const std::vector<double>& v) { intensities_ = v; }

    virtual std::string     getSolutionRepresentation() override;
    virtual emili::Solution* clone() override;
    /** Solo exige x >= 0; Dmin/Dmax son penalización blanda en el objetivo, no restricción dura. */
    virtual bool             isFeasible() override;
    virtual ~ImrtSolution() {}
};


/*---------------------------------------------------------------------------*
 *                         INITIAL SOLUTIONS                                 *
 *---------------------------------------------------------------------------*/

class ImrtInitialSolution : public emili::InitialSolution {
protected:
    ImrtProblem& problem_;
    virtual emili::Solution* generate() = 0;

public:
    explicit ImrtInitialSolution(ImrtProblem& p)
        : emili::InitialSolution(p), problem_(p) {}

    virtual emili::Solution* generateSolution()      override;
    virtual emili::Solution* generateEmptySolution() override;
};

/** All dimlets set to zero intensity. */
class ZeroInitialSolution : public ImrtInitialSolution {
protected:
    virtual emili::Solution* generate() override;
public:
    explicit ZeroInitialSolution(ImrtProblem& p) : ImrtInitialSolution(p) {}
};

/** All dimlets set to a fixed uniform intensity. */
class UniformInitialSolution : public ImrtInitialSolution {
protected:
    double intensity_;
    virtual emili::Solution* generate() override;
public:
    UniformInitialSolution(ImrtProblem& p, double intensity)
        : ImrtInitialSolution(p), intensity_(intensity) {}
};

/** Random intensities drawn from U[0, max_intensity]. */
class RandomInitialSolution : public ImrtInitialSolution {
protected:
    double max_intensity_;
    virtual emili::Solution* generate() override;
public:
    RandomInitialSolution(ImrtProblem& p, double max_intensity)
        : ImrtInitialSolution(p), max_intensity_(max_intensity) {}
};


/*---------------------------------------------------------------------------*
 *                           NEIGHBORHOODS                                   *
 *---------------------------------------------------------------------------*/

class ImrtNeighborhood : public emili::Neighborhood {
protected:
    ImrtProblem&     problem_;
    int              n_dimlets_;
    double           max_intensity_;
    std::vector<int> active_dimlets_; // global dimlet indices for active angles

    /** Populate active_dimlets_ from angle index list (empty = all). */
    void buildActiveDimlets(const std::vector<int>& angle_idxs);

public:
    explicit ImrtNeighborhood(ImrtProblem& p);
};

/**
 * SingleBeamletShift — changes one dimlet intensity by ±delta.
 * Neighbourhood size: 2 × n_dimlets.
 */
class SingleBeamletShift : public ImrtNeighborhood {
protected:
    double delta_;
    int    cur_dimlet_;
    int    cur_active_idx_; // position in active_dimlets_ (when active set is used)
    int    cur_dir_;
    bool   first_;
    int    last_dimlet_;
    double last_old_value_;

    virtual emili::Solution* computeStep(emili::Solution* step) override;
    virtual void reverseLastMove(emili::Solution* step) override;

public:
    SingleBeamletShift(ImrtProblem& p, double delta)
        : ImrtNeighborhood(p), delta_(delta)
        , cur_dimlet_(0), cur_active_idx_(0), cur_dir_(1), first_(true) {}

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int  size() override {
        return active_dimlets_.empty() ? 2 * n_dimlets_
                                       : 2 * (int)active_dimlets_.size();
    }
};

/**
 * BeamletSwap — swaps the intensities of two dimlets.
 * Neighbourhood size: n_dimlets × (n_dimlets − 1) / 2.
 * No es una operación neutra: cada dimlet tiene su propia columna dispersa
 * de influencia de dosis, así que mover la misma intensidad a otro dimlet
 * redistribuye la dosis sobre boxets distintos y cambia el objetivo.
 */
class BeamletSwap : public ImrtNeighborhood {
protected:
    int  i_, j_;           // global dimlet indices (full mode)
    int  i_act_, j_act_;   // indices into active_dimlets_ (active-angles mode)
    bool first_;
    int  last_i_, last_j_;

    virtual emili::Solution* computeStep(emili::Solution* step) override;
    virtual void reverseLastMove(emili::Solution* step) override;

public:
    explicit BeamletSwap(ImrtProblem& p)
        : ImrtNeighborhood(p), i_(0), j_(1), i_act_(0), j_act_(1), first_(true) {}

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int  size() override {
        if (active_dimlets_.empty())
            return n_dimlets_ * (n_dimlets_ - 1) / 2;
        int na = (int)active_dimlets_.size();
        return na * (na - 1) / 2;
    }
};


/*---------------------------------------------------------------------------*
 *                            PERTURBATION                                   *
 *---------------------------------------------------------------------------*/

/** Randomly reassigns k dimlet intensities to U[0, max_intensity]. */
class RandomBeamletPerturbation : public emili::Perturbation {
protected:
    ImrtProblem& problem_;
    int          k_;
    double       max_intensity_;

public:
    RandomBeamletPerturbation(ImrtProblem& p, int k, double max_intensity)
        : problem_(p), k_(k), max_intensity_(max_intensity) {}

    virtual emili::Solution* perturb(emili::Solution* solution) override;
};


/*---------------------------------------------------------------------------*
 *                             ACCEPTANCE                                    *
 *---------------------------------------------------------------------------*/

class ImrtImproveAccept : public emili::Acceptance {
public:
    virtual emili::Solution* accept(emili::Solution* current,
                                    emili::Solution* candidate) override;
};


/*---------------------------------------------------------------------------*
 *                            TERMINATION                                    *
 *---------------------------------------------------------------------------*/

class ImrtMaxIterations : public emili::Termination {
protected:
    int max_, cur_;
public:
    explicit ImrtMaxIterations(int max) : max_(max), cur_(0) {}
    virtual bool terminate(emili::Solution* c, emili::Solution* n) override;
    virtual void reset() override { cur_ = 0; }
};

class ImrtFeasibleTermination : public emili::Termination {
protected:
    double tolerance_;
public:
    explicit ImrtFeasibleTermination(double tol = 1e-6) : tolerance_(tol) {}
    virtual bool terminate(emili::Solution* c, emili::Solution* n) override;
    virtual void reset() override {}
};

/*---------------------------------------------------------------------------*
 *                           TABU MEMORY                                     *
 *---------------------------------------------------------------------------*/

/**
 * ImrtTabuMemory — circular buffer of intensity vectors.
 * Una solución es tabú si su vector de intensidades coincide exactamente
 * (operator==, sin tolerancia de punto flotante) con alguna entrada del
 * buffer — casi-duplicados de una perturbación no se detectan.
 */
class ImrtTabuMemory : public emili::TabuMemory {
protected:
    int tenure_;
    std::vector<std::vector<double>> memory_;
    int head_;
    int count_;

public:
    explicit ImrtTabuMemory(int tenure)
        : tenure_(tenure), memory_(tenure), head_(0), count_(0) {}

    virtual bool tabu_check(emili::Solution* solution) override;
    virtual void forbid(emili::Solution* solution) override;
    virtual void reset() override { head_ = 0; count_ = 0; }
};

} // namespace imrt
} // namespace emili

#endif // IMRT_H
