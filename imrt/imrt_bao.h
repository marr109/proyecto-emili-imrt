#ifndef IMRT_BAO_H
#define IMRT_BAO_H

#include "../emilibase.h"
#include "imrt_fmo.h"
#include "imrt_fmo_source.h"
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <deque>

namespace emili {
namespace imrt {

/*---------------------------------------------------------------------------*
 *                            BAO SOLUTION                                   *
 *                                                                           *
 * A BaoSolution represents a selection of K gantry angles (the outer BAO   *
 * decision) together with the FMO-optimal beamlet intensities for that     *
 * selection (the inner FMO result). NOTE: on this branch the FMO solver    *
 * is stubbed out — see imrt_fmo.h / ampl_gurobi/.                          *
 *---------------------------------------------------------------------------*/
class BaoSolution : public emili::Solution {
public:
    std::vector<int>    active_angles_; // K sorted angle indices (into inst.angles)
    std::vector<int>    angle_degrees_; // K actual degree values (for display only)
    std::vector<double> intensities_;   // n_beamlets optimal x* from FMO

    // Desglose por órgano del objetivo agregado (ver ImrtFmoSolver::solve /
    // FmoResult) — mismo orden que BaoProblem::ptvOrganNames()/oarOrganNames().
    // Vacíos hasta la primera evaluación (generateEmptySolution no evalúa).
    std::vector<double> ptv_underdose_sq_;
    std::vector<double> oar_overdose_sq_;

    BaoSolution(const std::vector<int>& angles, int n_beamlets)
        : emili::Solution(1e30)
        , active_angles_(angles)
        , intensities_(n_beamlets, 0.0)
    {}

    virtual const void* getRawData() const override { return this; }
    virtual void setRawData(const void* data) override {
        if (data == this) return;
        const BaoSolution* o = static_cast<const BaoSolution*>(data);
        active_angles_     = o->active_angles_;
        angle_degrees_     = o->angle_degrees_;
        intensities_       = o->intensities_;
        ptv_underdose_sq_  = o->ptv_underdose_sq_;
        oar_overdose_sq_   = o->oar_overdose_sq_;
    }

    virtual emili::Solution* clone() override;
    virtual std::string getSolutionRepresentation() override;
    virtual bool isFeasible() override { return true; }
    virtual ~BaoSolution() {}
};


/*---------------------------------------------------------------------------*
 *                             BAO PROBLEM                                   *
 *                                                                           *
 * BaoProblem wraps ImrtFmoSolver.  Evaluating a BaoSolution means calling   *
 * the FMO solver to optimise beamlet intensities for the given angle       *
 * subset, then storing the result back into the solution.                   *
 *---------------------------------------------------------------------------*/
class BaoProblem : public emili::Problem {
    struct CachedFmoResult {
        std::vector<double> intensities;
        double objective;
        std::vector<double> ptv_underdose_sq;
        std::vector<double> oar_overdose_sq;
    };

    std::unique_ptr<IFmoDataSource> source_;  // owned; ImrtFmoSolver holds a ref to this
    ImrtFmoSolver     fmo_;
    std::vector<int>  angle_degrees_;  // catalog: angle_degrees_[i] = real degree of candidate angle i
    int            K_;
    bool           verbose_;
    std::ofstream  csv_file_;
    int            eval_count_;
    std::map<std::vector<int>, CachedFmoResult> fmo_cache_;
    bool           last_eval_cached_;

public:
    BaoProblem(std::unique_ptr<IFmoDataSource> source, std::vector<int> angle_degrees, int K)
        : source_(std::move(source)), fmo_(*source_), angle_degrees_(std::move(angle_degrees))
        , K_(K), verbose_(false), eval_count_(0)
        , last_eval_cached_(false) {}

    virtual double calcObjectiveFunctionValue(emili::Solution& s) override;
    virtual double evaluateSolution(emili::Solution& s) override;
    virtual int    problemSize() override { return nAngles(); }

    const IFmoDataSource& getSource() const { return *source_; }
    int  K()                    const { return K_; }
    int  nAngles()               const { return (int)angle_degrees_.size(); }
    int  angleDegree(int idx)    const { return angle_degrees_[idx]; }
    int  nDimlets()               const { return source_->n_dimlets(); }

    // Nombres de órgano en el mismo orden que BaoSolution::ptv_underdose_sq_/
    // oar_overdose_sq_ — usados para armar el encabezado del CSV una sola vez.
    std::vector<std::string> ptvOrganNames() const { return fmo_.ptvOrganNames(); }
    std::vector<std::string> oarOrganNames() const { return fmo_.oarOrganNames(); }

    void setVerbose(bool v) { verbose_ = v; fmo_.setVerbose(v); }
    // Number of distinct angle sets solved by the FMO (real, non-cached solves).
    size_t realSolveCount()  const { return fmo_cache_.size(); }
    bool isReady()          const { return fmo_.isReady(); }
    void openCsvLog(const std::string& path);
};


/*---------------------------------------------------------------------------*
 *                        INITIAL SOLUTIONS                                  *
 *---------------------------------------------------------------------------*/

/** Select the first K angles (indices 0, 1, …, K-1). */
class FirstKAnglesInit : public emili::InitialSolution {
    BaoProblem& bao_;
    // true (default): restringe el catálogo a múltiplos de 5° — comportamiento
    // histórico post-2026-09-01, requerido por AngleShiftNeighborhood. false:
    // catálogo completo (0..nAngles()-1), sin filtro — para reproducir las
    // corridas "unrestricted" previas a esa fecha.
    bool restrict_to_5deg_grid_;
public:
    explicit FirstKAnglesInit(BaoProblem& p, bool restrict_to_5deg_grid = true)
        : emili::InitialSolution(p), bao_(p),
          restrict_to_5deg_grid_(restrict_to_5deg_grid) {}

    virtual emili::Solution* generateSolution()      override;
    virtual emili::Solution* generateEmptySolution() override;
};

/** Select K angles uniformly at random. */
class RandomKAnglesInit : public emili::InitialSolution {
    BaoProblem& bao_;
    // Ver comentario en FirstKAnglesInit: mismo significado y mismo default.
    bool restrict_to_5deg_grid_;
public:
    explicit RandomKAnglesInit(BaoProblem& p, bool restrict_to_5deg_grid = true)
        : emili::InitialSolution(p), bao_(p),
          restrict_to_5deg_grid_(restrict_to_5deg_grid) {}

    virtual emili::Solution* generateSolution()      override;
    virtual emili::Solution* generateEmptySolution() override;
};


/*---------------------------------------------------------------------------*
 *                          ANGLE-SWAP NEIGHBORHOOD                          *
 *                                                                           *
 * Remove one active angle and add one inactive angle.                       *
 * Size = K × (n_angles − K).                                               *
 *---------------------------------------------------------------------------*/
class AngleSwapNeighborhood : public emili::Neighborhood {
    BaoProblem& bao_;
    int         n_angles_;

    // State for iteration
    std::vector<int> base_angles_;     // active angles at begin()
    std::vector<int> inactive_list_;   // inactive angles at begin()
    int  cur_active_idx_;              // which active angle to remove
    int  cur_inactive_idx_;            // which inactive angle to add
    bool first_;

    virtual emili::Solution* computeStep(emili::Solution* step)  override;
    virtual void reverseLastMove(emili::Solution* step)           override;

public:
    explicit AngleSwapNeighborhood(BaoProblem& p)
        : bao_(p), n_angles_(p.nAngles())
        , cur_active_idx_(0), cur_inactive_idx_(0), first_(true)
    {}

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int size() override;
};

/*---------------------------------------------------------------------------*
 *                         ANGLE-SHIFT NEIGHBORHOOD                          *
 *                                                                           *
 * Mueve un ángulo activo ±`step` posiciones dentro del catálogo ordenado    *
 * por valor de grado real (no por índice crudo del array, que está en      *
 * orden lexicográfico de string). El desplazamiento es circular: después   *
 * del último ángulo (mayor grado) se vuelve al primero (0°), ya que el     *
 * gantry rota en un círculo continuo de 360°.                              *
 * Tamaño (cota superior) = 2 × K.                                          *
 *                                                                           *
 * Orden de recorrido de slots: por defecto (random_order_ = false) fijo,   *
 * 0, 1, ..., K-1 — sesga a First Improvement hacia el slot 0 (siempre se   *
 * revisa primero). Si random_order_ = true, begin() sortea una permutación *
 * de los slots (Fisher-Yates) y ese es el orden usado; el modo fijo sigue  *
 * siendo el default para no romper corridas/scripts existentes.            *
 *---------------------------------------------------------------------------*/
class AngleShiftNeighborhood : public emili::Neighborhood {
    BaoProblem& bao_;
    int         n_angles_;
    int         step_;
    bool        random_order_;       // false (default): recorre slots 0..K-1 en orden fijo.
                                      // true: begin() genera una permutación aleatoria de
                                      // slots (Fisher-Yates, mismo estilo que
                                      // AngleShiftMultiPerturbation::perturb) y recorre esa.
    bool        circular_order_;     // false (default): sin efecto sobre el orden fijo/random.
                                      // true: begin() arranca slot_order_ en el slot que cambió
                                      // en la ronda anterior (detectado por diff contra
                                      // base_angles_ previo) y rota circularmente desde ahí,
                                      // en vez de reiniciar siempre en el slot 0. Si
                                      // random_order_ y circular_order_ son ambos true,
                                      // random_order_ tiene prioridad (ya era el modo
                                      // establecido primero) y circular_order_ se ignora.

    // Permutación de índices de catálogo ordenados por grado ascendente,
    // y su lookup inverso (índice de catálogo -> posición en degree_order_)
    std::vector<int> degree_order_;
    std::vector<int> degree_rank_;

    // Estado para la iteración
    std::vector<int> base_angles_;   // ángulos activos al llamar begin()
    std::vector<int> slot_order_;    // orden de recorrido de slots: identidad (0..K-1) si
                                      // !random_order_, o permutación aleatoria si random_order_
    int  cur_active_idx_;            // posición actual dentro de slot_order_ (no índice crudo del slot)
    int  cur_dir_;                   // 0 = -step, 1 = +step
    bool first_;

    void buildDegreeOrder();

    virtual emili::Solution* computeStep(emili::Solution* step)  override;
    virtual void reverseLastMove(emili::Solution* step)           override;

public:
    explicit AngleShiftNeighborhood(BaoProblem& p, int step = 1, bool random_order = false,
                                     bool circular_order = false)
        : bao_(p), n_angles_(p.nAngles()), step_(step), random_order_(random_order)
        , circular_order_(circular_order)
        , cur_active_idx_(0), cur_dir_(0), first_(true)
    {}

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int size() override;
};

/*---------------------------------------------------------------------------*
 *                     ANGLE MULTI-SHIFT NEIGHBORHOOD                        *
 *                                                                           *
 * Versión "inclusiva" de AngleShiftNeighborhood: en vez de un único step,   *
 * recibe una lista de steps (p.ej. {5, 10}) y el vecindario de cada ángulo  *
 * activo acumula ±cada step de la lista — más vecinos, más caro de evaluar, *
 * pero ve más del catálogo por ronda. Con steps={s} es equivalente a        *
 * AngleShiftNeighborhood con step=s.                                        *
 * Tamaño (cota superior) = 2 × K × steps.size().                            *
 *---------------------------------------------------------------------------*/
class AngleMultiShiftNeighborhood : public emili::Neighborhood {
    BaoProblem& bao_;
    int         n_angles_;
    std::vector<int> steps_;

    // Permutación de índices de catálogo ordenados por grado ascendente,
    // y su lookup inverso (índice de catálogo -> posición en degree_order_)
    std::vector<int> degree_order_;
    std::vector<int> degree_rank_;

    // Estado para la iteración
    std::vector<int> base_angles_;   // ángulos activos al llamar begin()
    int  cur_active_idx_;            // qué ángulo activo se está desplazando
    int  cur_step_idx_;              // qué step de steps_ se está probando
    int  cur_dir_;                   // 0 = -step, 1 = +step
    bool first_;

    void buildDegreeOrder();

    virtual emili::Solution* computeStep(emili::Solution* step)  override;
    virtual void reverseLastMove(emili::Solution* step)           override;

public:
    explicit AngleMultiShiftNeighborhood(BaoProblem& p, std::vector<int> steps)
        : bao_(p), n_angles_(p.nAngles()), steps_(std::move(steps))
        , cur_active_idx_(0), cur_step_idx_(0), cur_dir_(0), first_(true)
    {}

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int size() override;
};

/*---------------------------------------------------------------------------*
 *                  RESTRICTED ANGLE-SWAP NEIGHBORHOODS                      *
 *                                                                           *
 * Cheaper subsets of AngleSwapNeighborhood (one active angle replaced by   *
 * one inactive catalog angle), materialised as an explicit move list      *
 * (slot, replacement catalog index) built at every begin():               *
 *   SAMPLE: min(n, K*(n_angles-K)) distinct swaps drawn uniformly without *
 *           replacement from the full swap neighbourhood, re-drawn at     *
 *           every begin() (Gutierrez & Cabrera-Guerrero 2018, N2).        *
 *   WINDOW: replacement within +-delta degrees (circular) of the slot's   *
 *           current angle; deterministic.                                  *
 *   FLIP:   replacement within +-delta degrees of (angle+180) mod 360;    *
 *           deterministic (Aleman et al. 2008, "flip").                    *
 * Replacements already active are never generated, so K is preserved and  *
 * no neighbour contains duplicated angles. Degrees that are not in the    *
 * catalog are skipped. By default (restrict_to_5deg_grid_ = true, same    *
 * convention and default as RandomKAnglesInit) only catalog angles that   *
 * are multiples of 5 degrees are candidates, so a run started on the 5-deg *
 * grid stays on it ("restricted" = 72 angles); with the optional builder  *
 * token "unrestricted" the full catalog is used.                          *
 * The slot is replaced in place (no sort), like                           *
 * AngleShiftNeighborhood; the FMO cache key is sorted in evaluateSolution. *
 * random() draws one move uniformly from the same definition (for SA).    *
 *---------------------------------------------------------------------------*/
class AngleSwapSubsetNeighborhood : public emili::Neighborhood {
public:
    enum Mode { SAMPLE, WINDOW, FLIP };

private:
    BaoProblem& bao_;
    int         n_angles_;
    Mode        mode_;
    int         param_;                // SAMPLE: n; WINDOW/FLIP: delta in degrees
    bool        restrict_to_5deg_grid_; // true: only catalog angles with degree % 5 == 0

    std::vector<int> degree_to_idx_;   // degree (0..359) -> catalog index, -1 if absent

    // State for iteration
    std::vector<int> base_angles_;               // active angles at begin()
    std::vector<std::pair<int,int>> moves_;      // (slot, replacement catalog index)
    int  cur_move_;
    bool first_;

    void buildDegreeLookup();
    bool isCandidate(int idx) const {
        return !restrict_to_5deg_grid_ || ((bao_.angleDegree(idx) % 360) + 360) % 5 == 0;
    }
    // All moves of this neighbourhood's definition for `angles` (SAMPLE: the
    // full swap neighbourhood, before sampling).
    std::vector<std::pair<int,int>> enumerateMoves(const std::vector<int>& angles);

    virtual emili::Solution* computeStep(emili::Solution* step)  override;
    virtual void reverseLastMove(emili::Solution* step)           override;

public:
    AngleSwapSubsetNeighborhood(BaoProblem& p, Mode mode, int param,
                                bool restrict_to_5deg_grid = true)
        : bao_(p), n_angles_(p.nAngles()), mode_(mode), param_(param)
        , restrict_to_5deg_grid_(restrict_to_5deg_grid)
        , cur_move_(0), first_(true)
    { buildDegreeLookup(); }

    virtual emili::Neighborhood::NeighborhoodIterator
            begin(emili::Solution* base) override;
    virtual void reset() override;
    virtual emili::Solution* random(emili::Solution* s) override;
    virtual int size() override;
};

/*---------------------------------------------------------------------------*
 *                  ANGLE SHIFT PERTURBATION, multi (for ILS)                *
 *                                                                           *
 * Perturbación descrita en la reunión con Leslie: desplaza numSteps         *
 * ángulos activos DISTINTOS (permutación sin reposición de los slots, no    *
 * numSteps sorteos independientes que pueden repetir el mismo slot), cada   *
 * uno una magnitud aleatoria en (step, 2*step) — nunca step exacto, para no *
 * quedar atrapado en la misma clase módulo step que usa el vecindario de    *
 * búsqueda local (ver AngleShiftNeighborhood::random).                      *
 *---------------------------------------------------------------------------*/
class AngleShiftMultiPerturbation : public emili::Perturbation {
    BaoProblem& bao_;
    int         n_angles_;
    int         step_;
    int         numSteps_;

    // Modo de magnitud discreta (opcional): si no esta vacio, cada angulo se
    // mueve +-m grados con m sorteado de esta lista (signo aleatorio) en vez
    // del rango (step, 2*step) posiciones. Vacio = modo rango (original).
    std::vector<int> magnitudes_;
    bool stepwise_;   // opt-in: evaluar tras cada movimiento individual
    bool replace_;    // opt-in: slots sorteados con reposición (un ángulo puede moverse más de una vez)

    std::vector<int> degree_order_;
    std::vector<int> degree_rank_;
    std::vector<int> degree_to_idx_;   // grado (0..359) -> indice de catalogo, -1 si no existe

    void buildDegreeOrder();

public:
    AngleShiftMultiPerturbation(BaoProblem& bao, int step, int numSteps,
                                const std::vector<int>& magnitudes = std::vector<int>(),
                                bool stepwise = false, bool replace = false)
        : emili::Perturbation(), bao_(bao), n_angles_(bao.nAngles())
        , step_(step), numSteps_(numSteps), magnitudes_(magnitudes), stepwise_(stepwise), replace_(replace) {}

    virtual emili::Solution* perturb(emili::Solution* solution) override;
};

/*---------------------------------------------------------------------------*
 *                 ADAPTIVE ANGLE SHIFT PERTURBATION (for ILS)               *
 *                                                                           *
 * Misma mecánica de movimiento que AngleShiftMultiPerturbation (permutación *
 * sin reposición de slots, magnitud aleatoria en (step,2*step) por slot),   *
 * pero progresiva: arranca moviendo 1 solo ángulo (nivel 0) y, si el        *
 * objetivo recibido en una llamada a perturb() no mejora respecto al        *
 * recibido en la llamada anterior durante stagnationThreshold llamadas      *
 * seguidas, sube un ángulo más por nivel (hasta maxNumSteps) — mismo        *
 * espíritu que MultiScaleAngleShake (k+1 intercambios para intensidad k),   *
 * aplicado acá al desplazamiento en vez del swap. La magnitud del           *
 * desplazamiento (step) se mantiene fija en todos los niveles; solo escala  *
 * la cantidad de ángulos movidos. En cuanto una llamada trae una mejora     *
 * real, vuelve a nivel 0 (1 ángulo). No hay ningún hook nuevo en el         *
 * framework genérico: perturb() ya recibe la solución actual aceptada por   *
 * el ILS en cada ronda, así que basta con comparar ese valor contra el      *
 * guardado en la llamada anterior (mismo patrón de auto-observación ya      *
 * usado en AngleShiftNeighborhood::begin() para el vecindario circular).    *
 *---------------------------------------------------------------------------*/
class AdaptiveAngleShiftPerturbation : public emili::Perturbation {
    BaoProblem& bao_;
    int         n_angles_;
    int         step_;
    int         maxNumSteps_;
    int         stagnationThreshold_;

    std::vector<int> magnitudes_;   // ver AngleShiftMultiPerturbation (vacio = modo rango)
    bool stepwise_;
    bool replace_;   // ver AngleShiftMultiPerturbation

    std::vector<int> degree_order_;
    std::vector<int> degree_rank_;
    std::vector<int> degree_to_idx_;

    double last_objective_;
    int    stagnation_count_;
    int    level_;   // numSteps efectivo = min(level_ + 1, maxNumSteps_)

    void buildDegreeOrder();

public:
    AdaptiveAngleShiftPerturbation(BaoProblem& bao, int step,
                                    int maxNumSteps, int stagnationThreshold,
                                    const std::vector<int>& magnitudes = std::vector<int>(),
                                    bool stepwise = false, bool replace = false)
        : emili::Perturbation(), bao_(bao), n_angles_(bao.nAngles())
        , step_(step), maxNumSteps_(maxNumSteps)
        , stagnationThreshold_(stagnationThreshold)
        , magnitudes_(magnitudes), stepwise_(stepwise), replace_(replace)
        , last_objective_(1e30), stagnation_count_(0), level_(0) {}

    virtual emili::Solution* perturb(emili::Solution* solution) override;
};

/*---------------------------------------------------------------------------*
 *                       ANGLE PERTURBATION (for ILS)                        *
 *                                                                           *
 * Randomly replaces p active angles with p inactive angles.                 *
 * Enables escaping local optima in an Iterated Local Search.                *
 *---------------------------------------------------------------------------*/
class RandomAnglesPerturbation : public emili::Perturbation {
    BaoProblem& bao_;
    int         p_;   // number of angles to swap out

public:
    RandomAnglesPerturbation(BaoProblem& bao, int p)
        : emili::Perturbation(), bao_(bao), p_(p) {}

    virtual emili::Solution* perturb(emili::Solution* current) override;
};


/*---------------------------------------------------------------------------*
 *                      BAO ACCEPTANCE (improve, ILS)                        *
 *                                                                           *
 * Standard "accept if strictly better" rule. When reject_repeated_ is set,  *
 * a candidate whose angle set was already seen earlier in the run is       *
 * rejected even if it improves on the current solution, so ILS keeps       *
 * diversifying instead of resettling on a known local optimum.             *
 *---------------------------------------------------------------------------*/
class BaoImproveAccept : public emili::Acceptance {
    bool                           reject_repeated_;
    std::set<std::vector<int>>     visited_;

public:
    explicit BaoImproveAccept(bool reject_repeated = false)
        : reject_repeated_(reject_repeated) {}

    virtual emili::Solution* accept(emili::Solution* current,
                                     emili::Solution* candidate) override;
    virtual void reset() override { visited_.clear(); }
};


/*---------------------------------------------------------------------------*
 *                    BAO TABU MEMORY (angle-set based)                       *
 *                                                                           *
 * Stores K-element active angle sets in a circular buffer.                  *
 * Suitable for the BAO level where moves are angle swaps.                   *
 *---------------------------------------------------------------------------*/
class BaoTabuMemory : public emili::TabuMemory {
protected:
    std::vector<std::vector<int>> memory_;
    int head_, count_;
    int tenure_;   // local copy (base field is private)

public:
    explicit BaoTabuMemory(int tenure)
        : emili::TabuMemory(tenure)
        , memory_(tenure)
        , head_(0), count_(0), tenure_(tenure) {}

    // Legacy semantic used by the generic EMILI `tabu` token: true when the
    // angle set IS in the tabu list. Kept unchanged for reproducibility.
    virtual bool tabu_check(emili::Solution* s) override { return isTabu(s); }
    virtual void forbid(emili::Solution* s)     override;
    virtual void reset() override { head_ = 0; count_ = 0; }

    // Unambiguous query: true if the angle set is stored in the tabu list.
    virtual bool isTabu(emili::Solution* s) const;

    int  getTenure() const { return tenure_; }

    // Angles currently banned by the memory (empty for configuration memories).
    virtual std::vector<int> bannedAngles() const { return {}; }
};


/*---------------------------------------------------------------------------*
 *                  ADAPTIVE BAO TABU MEMORY (reactive tenure)               *
 *                                                                           *
 * Extends BaoTabuMemory with dynamic tenure adjustment:                     *
 *   - Oscillation detected (revisit): tenure increases by 1 (up to max).   *
 *   - Quiet period (2×tenure steps without revisit): tenure decreases by 1. *
 *---------------------------------------------------------------------------*/
class AdaptiveBaoTabuMemory : public BaoTabuMemory {
    int tenure_min_, tenure_max_;
    int since_last_revisit_;

    void resizeTenure(int new_tenure);

public:
    AdaptiveBaoTabuMemory(int tenure_min, int tenure_max)
        : BaoTabuMemory((tenure_min + tenure_max) / 2)
        , tenure_min_(tenure_min), tenure_max_(tenure_max)
        , since_last_revisit_(0) {}

    virtual void forbid(emili::Solution* s) override;
};


/*---------------------------------------------------------------------------*
 *                  BAO ATTRIBUTE TABU MEMORY (token TBao_attr)               *
 *                                                                           *
 * Move-attribute memory instead of whole-solution memory. forbid() is       *
 * called with the new current solution after every move; the angles that    *
 * left the active set since the previous call are pushed into a FIFO of     *
 * `tenure` angles. A candidate is tabu iff it contains an angle that is in  *
 * the FIFO and not in the current solution, i.e. it would re-insert an      *
 * angle removed less than `tenure` removals ago. The comparison is on angle *
 * sets, so slot order is irrelevant.                                        *
 *---------------------------------------------------------------------------*/
class BaoAttributeTabuMemory : public BaoTabuMemory {
    std::deque<int> removed_;       // most recent at the back, at most tenure_ entries
    std::vector<int> current_;      // active angles of the last forbid() call
    bool has_current_;

public:
    explicit BaoAttributeTabuMemory(int tenure)
        : BaoTabuMemory(tenure), has_current_(false) {}

    virtual bool isTabu(emili::Solution* s) const override;
    virtual void forbid(emili::Solution* s) override;
    virtual void reset() override { BaoTabuMemory::reset(); removed_.clear(); current_.clear(); has_current_ = false; }
    virtual std::vector<int> bannedAngles() const override { return std::vector<int>(removed_.begin(), removed_.end()); }
};


/*---------------------------------------------------------------------------*
 *                     BAO TABU SEARCH (token `tabubao`)                     *
 *                                                                           *
 * Textbook tabu search, independent of emili::FirstTabuSearch /             *
 * BestTabuSearch (which rebuild the neighbourhood from bestSoFar):          *
 *   - neighbours are generated from the CURRENT solution;                   *
 *   - a neighbour is admissible if it is not tabu, or if it beats the best  *
 *     so far (aspiration);                                                  *
 *   - best:  move to the best admissible neighbour, even if worse;          *
 *     first: move to the first admissible neighbour improving the current   *
 *            one, otherwise to the best admissible neighbour of the scan;   *
 *   - if every neighbour is tabu, move to the best neighbour overall;       *
 *   - the visited solution is forbidden after every move;                   *
 *   - bestSoFar is updated separately (read by the -it SIGALRM handler).    *
 *   - the strategy token (first|best) is optional; the default is best.    *
 * Stops on the termination criterion or, if stall_limit > 0, after          *
 * stall_limit consecutive iterations without a real (non-cached) FMO solve. *
 *---------------------------------------------------------------------------*/
class BaoTabuSearch : public emili::LocalSearch {
    BaoTabuMemory& memory_;
    BaoProblem&    bao_;
    bool           first_;
    long           stall_limit_;
public:
    BaoTabuSearch(emili::InitialSolution& in, emili::Termination& te,
                  emili::Neighborhood& ne, BaoTabuMemory& mem, BaoProblem& bao,
                  bool first, long stall_limit)
        : emili::LocalSearch(in, te, ne), memory_(mem), bao_(bao)
        , first_(first), stall_limit_(stall_limit) {}

    virtual emili::Solution* search(emili::Solution* initial) override;
};


/*---------------------------------------------------------------------------*
 *                    MULTI-SCALE SHAKE (for rVNS)                           *
 *                                                                           *
 * shake(s, k) performs (k+1) random angle swaps.                           *
 * Provides shaking strengths k = 0, 1, …, p_max-1 for emili::GVNS.        *
 *---------------------------------------------------------------------------*/
class MultiScaleAngleShake : public emili::Shake {
    BaoProblem& bao_;

public:
    MultiScaleAngleShake(BaoProblem& bao, int p_max)
        : emili::Shake(p_max), bao_(bao) {}

    virtual emili::Solution* shake(emili::Solution* s, int k) override;
};


/*---------------------------------------------------------------------------*
 *                  GREEDY ANGLES PERTURBATION (Iterated Greedy)             *
 *                                                                           *
 * Destroys D randomly chosen active angles, then reconstructs them          *
 * greedily: at each step, every inactive candidate is evaluated via the    *
 * FMO solver and the one minimising the FMO objective is added.             *
 *---------------------------------------------------------------------------*/
class GreedyAnglesPerturbation : public emili::Perturbation {
    BaoProblem& bao_;
    int         D_;   // number of angles to destroy and rebuild

public:
    GreedyAnglesPerturbation(BaoProblem& bao, int D)
        : emili::Perturbation(), bao_(bao), D_(D) {}

    virtual emili::Solution* perturb(emili::Solution* current) override;
};

/*---------------------------------------------------------------------------*
 *              SAMPLED GREEDY ANGLES PERTURBATION (sampled IG)              *
 *                                                                           *
 * Like GreedyAnglesPerturbation, but each reconstruction step evaluates     *
 * only `sample` randomly drawn inactive angles (restricted to the 5-degree  *
 * grid when restrict5), so a perturbation costs at most D*sample FMO        *
 * solves (cache hits are free).                                             *
 *---------------------------------------------------------------------------*/
class SampledGreedyAnglesPerturbation : public emili::Perturbation {
    BaoProblem& bao_;
    int         D_;          // number of angles to destroy and rebuild
    int         sample_;     // candidates evaluated per reconstruction step
    bool        restrict5_;  // only 5-degree-grid angles are candidates

public:
    SampledGreedyAnglesPerturbation(BaoProblem& bao, int D, int sample, bool restrict5 = true)
        : emili::Perturbation(), bao_(bao), D_(D), sample_(sample), restrict5_(restrict5) {}

    virtual emili::Solution* perturb(emili::Solution* current) override;
};

} // namespace imrt
} // namespace emili

#endif // IMRT_BAO_H
