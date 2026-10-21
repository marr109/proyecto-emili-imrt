#include "imrt_builder.h"
#include "cerr_instance.h"
#include "imrt_fmo_source.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// ── Token keywords ────────────────────────────────────────────────────────────
#define PROBLEM_IMRT    "imrt"
#define PROBLEM_BAO     "baoimrt"
#define INIT_ZERO       "izero"
#define INIT_UNIFORM    "iuniform"
#define INIT_RANDOM     "irandom"
#define INIT_FIRSTK     "ifirstk"
#define INIT_RANDOMK    "irandomk"
#define INIT_K_UNRESTRICTED "unrestricted"
#define NEIGH_SHIFT     "nshift"
#define NEIGH_SWAP      "nswap"
#define NEIGH_ANGSWAP   "nangswap"
#define NEIGH_ANGSHIFT  "nangshift"
#define NEIGH_ANGSHIFT_MULTI "nangshiftmulti"
#define NEIGH_ANGSHIFT_RANDORDER "randorder"
#define NEIGH_ANGSHIFT_CIRCULARORDER "circular"
#define PERT_RANDOM     "prandom"
#define PERT_ANGSWAP    "prangswap"
#define PERT_ANGSHIFT   "prangshift"
#define PERT_ANGSHIFT_ADAPTIVE "prangshiftadaptive"
#define PERT_ANGSHIFT_MAG "prangshiftmag"
#define PERT_ANGSHIFT_ADAPTIVE_MAG "prangshiftadaptivemag"
#define ACC_IMPROVE     "aimprove"
#define ACC_BAO_IMPROVE "baoimprove"
#define ACC_REJECT_REPEATED "rejectrepeated"
#define TERM_MAXITER    "tmaxiter"
#define TERM_FEASIBLE   "tfeasible"
#define TERM_TIME       "ttime"
#define TABU_ALL_SOL    "Tabu_all_solution"
#define TABU_1_CAR      "Tabu_1_car"
#define TABU_1_POS      "Tabu_1_position"
#define OPT_NACTIVE        "nactive"
#define OPT_VERBOSE        "verbose"
#define OPT_CSV            "csv"
// ── New metaheuristic components ──────────────────────────────────────────────
#define PERT_GREEDY        "pgreedy"
#define TABU_BAO_FIXED     "TBao_fixed"
#define TABU_BAO_ADAPTIVE  "TBao_adaptive"
#define SHAKE_BANGSHAKE    "bangshake"

namespace prs {
namespace imrt {

/*---------------------------------------------------------------------------*
 * Helpers
 *---------------------------------------------------------------------------*/

emili::imrt::ImrtProblem* ImrtBuilder::castProblem()
{
    return static_cast<emili::imrt::ImrtProblem*>(gp.getInstance());
}

emili::imrt::BaoProblem* ImrtBuilder::castBaoProblem()
{
    return dynamic_cast<emili::imrt::BaoProblem*>(gp.getInstance());
}

bool ImrtBuilder::isBaoProblem()
{
    return castBaoProblem() != nullptr;
}

/*---------------------------------------------------------------------------*
 * Problem identification
 *---------------------------------------------------------------------------*/

bool ImrtBuilder::isCompatibleWith(char* problem_definition)
{
    if (strcmp(problem_definition, PROBLEM_IMRT) == 0) return true;
    if (strcmp(problem_definition, PROBLEM_BAO)  == 0) return true;
    return false;
}

bool ImrtBuilder::canOpenInstance(char* problem_definition)
{
    return isCompatibleWith(problem_definition);
}

/*---------------------------------------------------------------------------*
 * openInstance — dispatches on "imrt" vs "baoimrt"
 *
 * Token layout when called:
 *   tokens[0]            = executable path  (skipped by parseParams)
 *   tokens[1]            = instance directory
 *   tokens[currentToken] = problem keyword  ("imrt" or "baoimrt") ← peek here
 *   tokens[currentToken+1..] = optional K (BAO), then algorithm tokens
 *---------------------------------------------------------------------------*/

emili::Problem* ImrtBuilder::openInstance()
{
    // The problem keyword is the CURRENT token (not yet consumed).
    // peek() returns tokens[currentToken] without advancing.
    bool is_bao = false;
    int  K      = 0;

    {
        char* kw = tm.peek();
        if (kw && strcmp(kw, PROBLEM_BAO) == 0)
            is_bao = true;
    }

    // Instance directory lives at absolute token index 1 (argv[1]).
    char* dir = tm.tokenAt(1);
    tm.nextToken();   // consume the problem keyword ("imrt" or "baoimrt")
    prs::check(dir, "IMRT: expected path to instance directory");

    if (is_bao)
        K = tm.getInteger();   // consume K (number of gantry angles to select)

    // Only the raw CERR export layout (instances/CERR_Prostate) is supported
    // as an instance format on this branch — the CORT/old-format ImrtInstance
    // loader was removed once ImrtProblem stopped depending on it.
    if (!emili::imrt::CerrFmoSource::looksLikeCerrDir(dir)) {
        std::cerr << "[IMRT] " << dir
                  << " is not a CERR-format instance directory "
                     "(beamletIndex.txt not found).\n";
        exit(-1);
    }

    auto source = std::unique_ptr<emili::imrt::IFmoDataSource>(
        new emili::imrt::CerrFmoSource(dir));
    auto* cerr_src = static_cast<emili::imrt::CerrFmoSource*>(source.get());
    int n_angles = cerr_src->nAnglesTotal();

    // angle_idx IS the degree value: beamletIndex.txt covers a contiguous
    // 0..359 range at 1 degree spacing (verified against instances/CERR_Prostate).
    std::vector<int> angle_degrees(n_angles);
    for (int i = 0; i < n_angles; ++i) angle_degrees[i] = i;

    prs::printTabPlusOne("directory", dir);
    prs::printTab("CERR-format instance detected");
    prs::printTabPlusOne("angles",  n_angles);
    prs::printTabPlusOne("dimlets", source->n_dimlets());

    if (is_bao) {
        if (K <= 0 || K > n_angles) {
            std::cerr << "[BAO] K=" << K << " invalido para "
                      << n_angles << " angulos\n";
            exit(-1);
        }
        prs::printTab("BAO problem (angle search)");
        prs::printTabPlusOne("K (angulos activos)", K);

        emili::imrt::BaoProblem* prob =
            new emili::imrt::BaoProblem(std::move(source), angle_degrees, K);
        if (!prob->isReady()) {
            std::cerr << "[BAO] FMO solver initialization failed\n";
            exit(-1);
        }

        if (tm.checkToken(OPT_VERBOSE)) {
            prob->setVerbose(true);
            prs::printTab("modo verbose BAO activado");
        }
        if (tm.checkToken(OPT_CSV)) {
            char* path = tm.nextToken();
            if (path) {
                prob->openCsvLog(std::string(path));
                prs::printTabPlusOne("CSV log", path);
            }
        }
        return prob;
    }

    // Classic FMO problem
    prs::printTab("IMRT problem loaded");
    emili::imrt::ImrtProblem* prob =
        new emili::imrt::ImrtProblem(std::move(source), angle_degrees);

    if (tm.checkToken(OPT_NACTIVE)) {
        int k = tm.getInteger();
        prob->setActiveAngles(k);
        prs::printTabPlusOne("angulos activos", k);
    }
    if (tm.checkToken(OPT_VERBOSE)) {
        prob->setVerboseFmo(true);
        prs::printTab("modo verbose FMO activado");
    }
    return prob;
}

/*---------------------------------------------------------------------------*
 * buildInitialSolution
 *---------------------------------------------------------------------------*/

emili::InitialSolution* ImrtBuilder::buildInitialSolution()
{
    prs::incrementTabLevel();
    emili::InitialSolution* init = nullptr;

    if (isBaoProblem()) {
        emili::imrt::BaoProblem* prob = castBaoProblem();
        if (tm.checkToken(INIT_FIRSTK)) {
            // Token opcional final "unrestricted": desactiva el filtro a
            // múltiplos de 5° y usa el catálogo completo de ángulos. Sin
            // este token, comportamiento idéntico al histórico: catálogo
            // restringido a la reja de 5° — ningún script existente cambia
            // de resultado.
            bool unrestricted = tm.checkToken(INIT_K_UNRESTRICTED);
            prs::printTab("BAO initial solution: first K angles");
            prs::printTabPlusOne("catalog", unrestricted ? "unrestricted (full 360)" : "restricted (mod 5)");
            init = new emili::imrt::FirstKAnglesInit(*prob, !unrestricted);
        }
        else if (tm.checkToken(INIT_RANDOMK)) {
            // Ver comentario en INIT_FIRSTK: mismo token, mismo significado.
            bool unrestricted = tm.checkToken(INIT_K_UNRESTRICTED);
            prs::printTab("BAO initial solution: random K angles");
            prs::printTabPlusOne("catalog", unrestricted ? "unrestricted (full 360)" : "restricted (mod 5)");
            init = new emili::imrt::RandomKAnglesInit(*prob, !unrestricted);
        }
        prs::decrementTabLevel();
        return init;
    }

    emili::imrt::ImrtProblem* prob = castProblem();

    if (tm.checkToken(INIT_ZERO)) {
        prs::printTab("initial solution: zero");
        init = new emili::imrt::ZeroInitialSolution(*prob);
    }
    else if (tm.checkToken(INIT_UNIFORM)) {
        double x = tm.getDecimal();
        prs::printTab("initial solution: uniform");
        prs::printTabPlusOne("intensity", x);
        init = new emili::imrt::UniformInitialSolution(*prob, x);
    }
    else if (tm.checkToken(INIT_RANDOM)) {
        double max = tm.getDecimal();
        prs::printTab("initial solution: random");
        prs::printTabPlusOne("max_intensity", max);
        init = new emili::imrt::RandomInitialSolution(*prob, max);
    }

    prs::decrementTabLevel();
    return init;
}

/*---------------------------------------------------------------------------*
 * buildNeighborhood
 *---------------------------------------------------------------------------*/

emili::Neighborhood* ImrtBuilder::buildNeighborhood()
{
    prs::incrementTabLevel();
    emili::Neighborhood* neigh = nullptr;

    if (isBaoProblem()) {
        emili::imrt::BaoProblem* prob = castBaoProblem();
        if (tm.checkToken(NEIGH_ANGSWAP)) {
            prs::printTab("BAO neighborhood: angle swap");
            neigh = new emili::imrt::AngleSwapNeighborhood(*prob);
        }
        else if (tm.checkToken(NEIGH_ANGSHIFT)) {
            int step = tm.getInteger();
            // Token opcional final "randorder": habilita el orden de slots
            // aleatorio (Fisher-Yates en begin(), ver AngleShiftNeighborhood).
            // Token opcional final alternativo "circular": arranca el recorrido
            // de slots en el slot recién modificado por la ronda anterior y
            // rota circularmente desde ahí (ver AngleShiftNeighborhood::begin()).
            // Sin ninguno de estos tokens, comportamiento idéntico al histórico:
            // orden fijo 0..K-1 — ningún script existente cambia de resultado.
            bool random_order = tm.checkToken(NEIGH_ANGSHIFT_RANDORDER);
            bool circular_order = tm.checkToken(NEIGH_ANGSHIFT_CIRCULARORDER);
            prs::printTab("BAO neighborhood: angle shift");
            prs::printTabPlusOne("step", step);
            prs::printTabPlusOne("slot order",
                random_order ? "random" : (circular_order ? "circular" : "fixed"));
            neigh = new emili::imrt::AngleShiftNeighborhood(*prob, step, random_order, circular_order);
        }
        else if (tm.checkToken(NEIGH_ANGSHIFT_MULTI)) {
            // Sintaxis: nangshiftmulti <n_steps> <step_1> ... <step_n>
            // Vecindario "inclusivo": acumula ±step_i para cada step_i de la lista,
            // en vez de un único step (ver AngleMultiShiftNeighborhood).
            int n_steps = tm.getInteger();
            std::vector<int> steps;
            steps.reserve(n_steps);
            for (int i = 0; i < n_steps; ++i) steps.push_back(tm.getInteger());
            prs::printTab("BAO neighborhood: angle shift (multi, inclusive)");
            for (int s : steps) prs::printTabPlusOne("step", s);
            neigh = new emili::imrt::AngleMultiShiftNeighborhood(*prob, steps);
        }
        prs::decrementTabLevel();
        return neigh;
    }

    emili::imrt::ImrtProblem* prob = castProblem();

    if (tm.checkToken(NEIGH_SHIFT)) {
        double delta = tm.getDecimal();
        prs::printTab("neighborhood: single beamlet shift");
        prs::printTabPlusOne("delta", delta);
        neigh = new emili::imrt::SingleBeamletShift(*prob, delta);
    }
    else if (tm.checkToken(NEIGH_SWAP)) {
        prs::printTab("neighborhood: beamlet swap");
        neigh = new emili::imrt::BeamletSwap(*prob);
    }

    prs::decrementTabLevel();
    return neigh;
}

// Reads "<count> <m1> ... <m_count>" for the discrete-magnitude perturbations.
// Fails loudly (exit) on count < 1 or any magnitude < 1.
static std::vector<int> readMagnitudes(TokenManager& tm)
{
    int count = tm.getInteger();
    if (count < 1) {
        std::cerr << "[IMRT] discrete-magnitude perturbation: count must be >= 1 (got "
                  << count << ")" << std::endl;
        exit(-1);
    }
    std::vector<int> mags;
    for (int i = 0; i < count; ++i) {
        int m = tm.getInteger();
        if (m < 1) {
            std::cerr << "[IMRT] discrete-magnitude perturbation: magnitude must be >= 1 (got "
                      << m << ")" << std::endl;
            exit(-1);
        }
        mags.push_back(m);
    }
    return mags;
}

static void printStepwise(bool sw)
{
    if (sw) prs::printTabPlusOne("stepwise", "evaluate after each individual move");
}

static void printMagnitudes(const std::vector<int>& mags)
{
    std::string s;
    for (size_t i = 0; i < mags.size(); ++i) s += (i ? "," : "") + std::to_string(mags[i]);
    prs::printTabPlusOne("magnitudes (deg)", s.c_str());
}

/*---------------------------------------------------------------------------*
 * buildPerturbation
 *---------------------------------------------------------------------------*/

emili::Perturbation* ImrtBuilder::buildPerturbation()
{
    prs::incrementTabLevel();
    emili::Perturbation* pert = nullptr;

    if (isBaoProblem()) {
        emili::imrt::BaoProblem* prob = castBaoProblem();
        if (tm.checkToken(PERT_ANGSWAP)) {
            int p = tm.getInteger();
            prs::printTab("BAO perturbation: random angle swap");
            prs::printTabPlusOne("swaps", p);
            pert = new emili::imrt::RandomAnglesPerturbation(*prob, p);
        }
        else if (tm.checkToken(PERT_GREEDY)) {
            int D = tm.getInteger();
            prs::printTab("BAO perturbation: greedy construction (Iterated Greedy)");
            prs::printTabPlusOne("D (destroy/rebuild)", D);
            pert = new emili::imrt::GreedyAnglesPerturbation(*prob, D);
        }
        else if (tm.checkToken(PERT_ANGSHIFT)) {
            // Sintaxis: prangshift <step> <numSteps>
            // Perturbación descrita por Leslie: desplaza numSteps ángulos
            // activos DISTINTOS (no numSteps sorteos con reposición, que
            // pueden repetir el mismo slot y terminar moviendo un solo
            // ángulo), cada uno una magnitud aleatoria en (step, 2*step) —
            // ver AngleShiftMultiPerturbation.
            int step      = tm.getInteger();
            int numSteps  = tm.getInteger();
            prs::printTab("BAO perturbation: multi angle shift (distinct slots)");
            prs::printTabPlusOne("step", step);
            prs::printTabPlusOne("numSteps", numSteps);
            bool sw = tm.checkToken("stepwise");
            printStepwise(sw);
            pert = new emili::imrt::AngleShiftMultiPerturbation(*prob, step, numSteps,
                                                                std::vector<int>(), sw);
        }
        else if (tm.checkToken(PERT_ANGSHIFT_ADAPTIVE)) {
            // Sintaxis: prangshiftadaptive <step> <maxNumSteps> <stagnationThreshold>
            // Misma mecánica de movimiento que AngleShiftMultiPerturbation,
            // pero progresiva: arranca moviendo 1 ángulo (nivel 0) y sube un
            // ángulo más por nivel cuando el objetivo no mejora entre
            // llamadas sucesivas a perturb() durante stagnationThreshold
            // rondas seguidas (hasta maxNumSteps); vuelve a nivel 0 apenas
            // hay una mejora real. La magnitud (step) queda fija en todos
            // los niveles. Ver AdaptiveAngleShiftPerturbation.
            int step                 = tm.getInteger();
            int maxNumSteps          = tm.getInteger();
            int stagnationThreshold  = tm.getInteger();
            prs::printTab("BAO perturbation: adaptive multi angle shift");
            prs::printTabPlusOne("step", step);
            prs::printTabPlusOne("maxNumSteps", maxNumSteps);
            prs::printTabPlusOne("stagnationThreshold", stagnationThreshold);
            bool sw = tm.checkToken("stepwise");
            printStepwise(sw);
            pert = new emili::imrt::AdaptiveAngleShiftPerturbation(
                *prob, step, maxNumSteps, stagnationThreshold, std::vector<int>(), sw);
        }
        else if (tm.checkToken(PERT_ANGSHIFT_MAG)) {
            // Sintaxis: prangshiftmag <numSteps> <count> <m1> ... <m_count>
            // Igual que prangshift pero cada angulo se mueve +-m grados, con m
            // sorteado de la lista y signo aleatorio (step no se usa).
            int numSteps = tm.getInteger();
            std::vector<int> mags = readMagnitudes(tm);
            prs::printTab("BAO perturbation: multi angle shift, discrete magnitudes");
            prs::printTabPlusOne("numSteps", numSteps);
            printMagnitudes(mags);
            bool sw = tm.checkToken("stepwise");
            printStepwise(sw);
            pert = new emili::imrt::AngleShiftMultiPerturbation(*prob, 1, numSteps, mags, sw);
        }
        else if (tm.checkToken(PERT_ANGSHIFT_ADAPTIVE_MAG)) {
            // Sintaxis: prangshiftadaptivemag <maxNumSteps> <stagnationThreshold> <count> <m1> ... <m_count>
            int maxNumSteps          = tm.getInteger();
            int stagnationThreshold  = tm.getInteger();
            std::vector<int> mags = readMagnitudes(tm);
            prs::printTab("BAO perturbation: adaptive multi angle shift, discrete magnitudes");
            prs::printTabPlusOne("maxNumSteps", maxNumSteps);
            prs::printTabPlusOne("stagnationThreshold", stagnationThreshold);
            printMagnitudes(mags);
            bool sw = tm.checkToken("stepwise");
            printStepwise(sw);
            pert = new emili::imrt::AdaptiveAngleShiftPerturbation(
                *prob, 1, maxNumSteps, stagnationThreshold, mags, sw);
        }
        prs::decrementTabLevel();
        return pert;
    }

    emili::imrt::ImrtProblem* prob = castProblem();

    if (tm.checkToken(PERT_RANDOM)) {
        int    k   = tm.getInteger();
        double max = tm.getDecimal();
        prs::printTab("perturbation: random beamlet");
        prs::printTabPlusOne("k", k);
        prs::printTabPlusOne("max_intensity", max);
        pert = new emili::imrt::RandomBeamletPerturbation(*prob, k, max);
    }

    prs::decrementTabLevel();
    return pert;
}

/*---------------------------------------------------------------------------*
 * buildAcceptance
 *---------------------------------------------------------------------------*/

emili::Acceptance* ImrtBuilder::buildAcceptance()
{
    prs::incrementTabLevel();
    emili::Acceptance* acc = nullptr;

    if (tm.checkToken(ACC_IMPROVE)) {
        prs::printTab("acceptance: improve");
        acc = new emili::imrt::ImrtImproveAccept();
    }
    else if (isBaoProblem() && tm.checkToken(ACC_BAO_IMPROVE)) {
        // Sintaxis: baoimprove [rejectrepeated]
        // Igual regla "improve" de siempre, pero consciente del conjunto de
        // ángulos activos (BaoSolution). rejectrepeated es opcional: si está
        // presente, rechaza un candidato que mejore pero repita un conjunto
        // de ángulos ya visitado en la corrida — ver BaoImproveAccept.
        bool reject_repeated = tm.checkToken(ACC_REJECT_REPEATED);
        prs::printTab("BAO acceptance: improve");
        prs::printTabPlusOne("reject repeated solutions", reject_repeated ? "true" : "false");
        acc = new emili::imrt::BaoImproveAccept(reject_repeated);
    }

    prs::decrementTabLevel();
    return acc;
}

/*---------------------------------------------------------------------------*
 * buildTermination
 *---------------------------------------------------------------------------*/

emili::Termination* ImrtBuilder::buildTermination()
{
    prs::incrementTabLevel();
    emili::Termination* term = nullptr;

    if (tm.checkToken(TERM_MAXITER)) {
        int n = tm.getInteger();
        prs::printTab("termination: max iterations");
        prs::printTabPlusOne("max", n);
        term = new emili::imrt::ImrtMaxIterations(n);
    }
    else if (tm.checkToken(TERM_FEASIBLE)) {
        prs::printTab("termination: feasible solution found");
        term = new emili::imrt::ImrtFeasibleTermination();
    }
    else if (tm.checkToken(TERM_TIME)) {
        // Reutiliza emili::TimedTermination (ya existente en el framework
        // genérico, wireada en generalParser.cpp bajo el token "time" para
        // el resto de EMILI) -- mide segundos de reloj desde reset() hasta
        // que terminate() detecta que se cumplió el máximo.
        float secs = (float)tm.getDecimal();
        prs::printTab("termination: max time");
        prs::printTabPlusOne("max_secs", secs);
        term = new emili::TimedTermination(secs);
    }

    prs::decrementTabLevel();
    return term;
}

/*---------------------------------------------------------------------------*
 * buildTabuTenure
 *---------------------------------------------------------------------------*/

emili::TabuMemory* ImrtBuilder::buildTabuTenure()
{
    prs::incrementTabLevel();
    emili::TabuMemory* mem = nullptr;

    if (tm.checkToken(TABU_ALL_SOL) ||
        tm.checkToken(TABU_1_CAR)   ||
        tm.checkToken(TABU_1_POS))
    {
        int tenure = tm.getInteger();
        prs::printTab("tabu memory: solution vector (circular buffer)");
        prs::printTabPlusOne("tenure", tenure);
        mem = new emili::imrt::ImrtTabuMemory(tenure);
    }
    else if (tm.checkToken(TABU_BAO_FIXED)) {
        int tenure = tm.getInteger();
        prs::printTab("BAO tabu memory: fixed tenure on angle sets");
        prs::printTabPlusOne("tenure", tenure);
        mem = new emili::imrt::BaoTabuMemory(tenure);
    }
    else if (tm.checkToken(TABU_BAO_ADAPTIVE)) {
        int tmin = tm.getInteger();
        int tmax = tm.getInteger();
        prs::printTab("BAO tabu memory: adaptive tenure on angle sets");
        prs::printTabPlusOne("tenure_min", tmin);
        prs::printTabPlusOne("tenure_max", tmax);
        mem = new emili::imrt::AdaptiveBaoTabuMemory(tmin, tmax);
    }

    prs::decrementTabLevel();
    return mem;
}

/*---------------------------------------------------------------------------*
 * buildShake — rVNS multi-scale shake for BAO
 *---------------------------------------------------------------------------*/

emili::Shake* ImrtBuilder::buildShake()
{
    prs::incrementTabLevel();
    emili::Shake* sh = nullptr;

    if (isBaoProblem()) {
        emili::imrt::BaoProblem* prob = castBaoProblem();
        if (tm.checkToken(SHAKE_BANGSHAKE)) {
            int p_max = tm.getInteger();
            prs::printTab("BAO shake: multi-scale angle swap (rVNS)");
            prs::printTabPlusOne("p_max", p_max);
            sh = new emili::imrt::MultiScaleAngleShake(*prob, p_max);
        }
        prs::decrementTabLevel();
        return sh;
    }

    prs::decrementTabLevel();
    return sh;
}

} // namespace imrt
} // namespace prs
