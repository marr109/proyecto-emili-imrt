#include "imrt_fmo.h"

#include <ampl/ampl.h>

#include <chrono>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef EMILI_REPO_ROOT
#define EMILI_REPO_ROOT "."
#endif

namespace emili {
namespace imrt {

namespace {

std::string envOr(const char* var, const std::string& fallback)
{
    const char* v = std::getenv(var);
    return (v && *v) ? std::string(v) : fallback;
}

std::string repoPath(const std::string& rel)
{
    return std::string(EMILI_REPO_ROOT) + "/" + rel;
}

} // namespace

/*---------------------------------------------------------------------------*
 * Constructor / destructor
 *---------------------------------------------------------------------------*/
ImrtFmoSolver::ImrtFmoSolver(const IFmoDataSource& source)
    : source_(source), ready_(false), n_ptv_(0), n_oar_(0)
{
    precompute();
    try {
        initAmpl();
        ready_ = true;
    } catch (const std::exception& e) {
        std::cerr << "[FMO] ERROR: AMPL/Gurobi initialization failed: " << e.what() << "\n";
        ready_ = false;
    }
}

// Out-of-line so ~unique_ptr<ampl::AMPL> only needs the complete type here,
// not wherever ImrtFmoSolver is used as a member (imrt_bao.h stays free of
// the AMPL headers).
ImrtFmoSolver::~ImrtFmoSolver() = default;

/*---------------------------------------------------------------------------*
 * precompute — organ boundaries (constant across solves for a given source)
 *---------------------------------------------------------------------------*/
void ImrtFmoSolver::precompute()
{
    n_ptv_ = 0;
    n_oar_ = 0;
    dmin_.clear();
    dmax_.clear();
    dmax_ptv_.clear();
    ptv_bounds_.clear();
    oar_bounds_.clear();

    for (const FmoOrganRef& o : source_.ptvOrgans()) {
        ptv_bounds_.push_back({o.name, n_ptv_, o.n_boxets});
        n_ptv_ += o.n_boxets;
        for (int b = 0; b < o.n_boxets; ++b) {
            dmin_.push_back(o.dmin);
            dmax_ptv_.push_back(o.dmax_ptv);
        }
    }
    for (const FmoOrganRef& o : source_.oarOrgans()) {
        oar_bounds_.push_back({o.name, n_oar_, o.n_boxets});
        n_oar_ += o.n_boxets;
        for (int b = 0; b < o.n_boxets; ++b)
            dmax_.push_back(o.dmax);
    }
}

/*---------------------------------------------------------------------------*
 * printSolveDims — advisor-requested sanity check (see git history / project
 * notes): print each organ's boxet count (fixed) and how many sparse dose
 * entries feed this specific solve (varies with the active-angle set), so a
 * silently-stale or repeated matrix is visible immediately instead of
 * hiding inside an otherwise-plausible objective value.
 *---------------------------------------------------------------------------*/
void ImrtFmoSolver::printSolveDims(int n_active_beamlets,
                                    const std::vector<int>& ptv_nnz,
                                    const std::vector<int>& oar_nnz) const
{
    std::cout << "  [FMO] beamlets_activos=" << n_active_beamlets << "\n";
    for (size_t i = 0; i < ptv_bounds_.size(); ++i)
        std::cout << "    PTV " << ptv_bounds_[i].name
                   << ": boxets=" << ptv_bounds_[i].n_boxets
                   << " entradas_dosis=" << ptv_nnz[i] << "\n";
    for (size_t i = 0; i < oar_bounds_.size(); ++i)
        std::cout << "    OAR " << oar_bounds_[i].name
                   << ": boxets=" << oar_bounds_[i].n_boxets
                   << " entradas_dosis=" << oar_nnz[i] << "\n";
    // Gurobi runs as a child process writing straight to the inherited
    // stdout fd; without an explicit flush here, this block-buffered
    // std::cout output can appear AFTER that child's output in the
    // terminal even though it was written first -- exactly backwards from
    // what the "print before solving" sanity check is meant to show.
    std::cout.flush();
}

/*---------------------------------------------------------------------------*
 * initAmpl — one-time environment/model setup + the static (per-source,
 * not per-solve) parameters.
 *---------------------------------------------------------------------------*/
void ImrtFmoSolver::initAmpl()
{
    std::string bin_dir = envOr("EMILI_AMPL_BIN_DIR",
        repoPath("ampl_gurobi/.venv/lib/python3.9/site-packages/ampl_module_base/bin"));
    std::string gurobi_bin = envOr("EMILI_GUROBI_BIN",
        repoPath("ampl_gurobi/.venv/lib/python3.9/site-packages/ampl_module_gurobi/bin/gurobi"));

    ampl::Environment env(bin_dir);
    ampl_.reset(new ampl::AMPL(env));
    ampl_->setOption("solver", gurobi_bin);

    // Solver options via EMILI_GUROBI_OPTIONS, passed through to gurobi_options
    // (e.g. "threads=8 barconvtol=1e-6 outlev=1"). Setting the shell environment
    // variable does not reach the solver: ampl::Environment gives the AMPL child
    // its own environment, so the option has to be set on the session.
    std::string gurobi_opts = envOr("EMILI_GUROBI_OPTIONS", "");
    if (!gurobi_opts.empty()) {
        ampl_->setOption("gurobi_options", gurobi_opts);
        std::cout << "  gurobi_options: " << gurobi_opts << "\n";
    }

    // Model selection: EMILI_FMO_MODEL unset/"baseline" -> fmo.mod (single
    // w_over); "organs" -> fmo_organs.mod (one overdose weight per OAR).
    const std::string model_name = envOr("EMILI_FMO_MODEL", "baseline");
    if (model_name != "baseline" && model_name != "organs")
        throw std::runtime_error("EMILI_FMO_MODEL must be 'baseline' or 'organs', got '"
                                 + model_name + "'");
    const bool per_organ = (model_name == "organs");
    if (!per_organ) std::cout << "  fmo_model: baseline (w_over=" << source_.w_over() << ")\n";
    ampl_->read(repoPath(per_organ ? "ampl_gurobi/fmo_organs.mod"
                                   : "ampl_gurobi/fmo.mod"));

    ampl_->getParameter("n_ptv").set(n_ptv_);
    ampl_->getParameter("n_oar").set(n_oar_);
    ampl_->getParameter("max_intensity").set(source_.max_intensity());
    ampl_->getParameter("w_under").set(source_.w_under());
    if (!per_organ) {
        ampl_->getParameter("w_over").set(source_.w_over());
    } else {
        const size_t n_org = oar_bounds_.size();
        if (n_org == 0)
            throw std::runtime_error("EMILI_FMO_MODEL=organs requires at least one OAR");

        // Weights: comma-separated, one per organ in oar_bounds_ order.
        std::vector<double> w_organ;
        const std::string w_str = envOr("EMILI_W_OVER_ORGANS", "");
        if (w_str.empty()) {
            w_organ.assign(n_org, source_.w_over());
        } else {
            std::stringstream ss(w_str);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                try {
                    size_t pos = 0;
                    w_organ.push_back(std::stod(tok, &pos));
                    if (tok.find_first_not_of(" \t", pos) != std::string::npos)
                        throw std::invalid_argument(tok);
                } catch (const std::exception&) {
                    throw std::runtime_error("EMILI_W_OVER_ORGANS: invalid number '" + tok + "'");
                }
            }
        }
        if (w_organ.size() != n_org)
            throw std::runtime_error("EMILI_W_OVER_ORGANS has " + std::to_string(w_organ.size())
                + " values but there are " + std::to_string(n_org) + " OAR organs");

        std::vector<double> oar_organ(n_oar_, 0.0);
        for (size_t oi = 0; oi < n_org; ++oi)
            for (int r = oar_bounds_[oi].row_off;
                 r < oar_bounds_[oi].row_off + oar_bounds_[oi].n_boxets; ++r)
                oar_organ[r] = (double)oi;

        const bool normalize = (envOr("EMILI_FMO_NORMALIZE", "0") == "1");
        ampl_->getParameter("n_organs").set((int)n_org);
        ampl_->getParameter("w_over_organ").setValues(ampl::Args(w_organ.data()), w_organ.size());
        if (n_oar_ > 0)
            ampl_->getParameter("oar_organ").setValues(ampl::Args(oar_organ.data()), oar_organ.size());
        ampl_->getParameter("normalize").set(normalize ? 1 : 0);

        std::cout << "  fmo_model: organs (normalize=" << (normalize ? 1 : 0) << ")\n";
        for (size_t oi = 0; oi < n_org; ++oi)
            std::cout << "    w_over_organ[" << oar_bounds_[oi].name << "]=" << w_organ[oi] << "\n";
        std::cout.flush();
    }
    ampl_->getParameter("w_ptv_over").set(source_.w_ptv_over());

    // dmin/dmax/dmax_ptv are indexed by PTV_B/OAR_B := 0..n-1, the same
    // ascending row order dmin_/dmax_/dmax_ptv_ were built in, so a plain
    // positional setValues matches AMPL's iteration order for these ranges.
    if (n_ptv_ > 0) {
        ampl_->getParameter("dmin").setValues(ampl::Args(dmin_.data()), dmin_.size());
        ampl_->getParameter("dmax_ptv").setValues(ampl::Args(dmax_ptv_.data()), dmax_ptv_.size());
    }
    if (n_oar_ > 0) {
        ampl_->getParameter("dmax").setValues(ampl::Args(dmax_.data()), dmax_.size());
    }
}

/*---------------------------------------------------------------------------*
 * Nombres de órgano, mismo orden que precompute() usó para construir
 * ptv_bounds_/oar_bounds_ (y por lo tanto el mismo orden en que solve()
 * agrega ptv_underdose_sq/oar_overdose_sq).
 *---------------------------------------------------------------------------*/
std::vector<std::string> ImrtFmoSolver::ptvOrganNames() const
{
    std::vector<std::string> names;
    names.reserve(ptv_bounds_.size());
    for (const auto& b : ptv_bounds_) names.push_back(b.name);
    return names;
}

std::vector<std::string> ImrtFmoSolver::oarOrganNames() const
{
    std::vector<std::string> names;
    names.reserve(oar_bounds_.size());
    for (const auto& b : oar_bounds_) names.push_back(b.name);
    return names;
}

/*---------------------------------------------------------------------------*
 * solve — refresh the per-solve data (active dimlets, sparse dose sets) and
 * resolve with Gurobi.
 *---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*
 * AlarmGuard — bloquea SIGALRM mientras dura el solve.
 *
 * El presupuesto de tiempo (-it) se entrega por SIGALRM, y su handler
 * (emili::finalise) imprime el reporte clinico: reserva memoria y usa
 * iostreams, nada de eso es async-signal-safe. El 95% del tiempo de reloj se
 * va dentro de este solve, casi todo en asignaciones para las ~850k tuplas de
 * dosis, asi que la señal caia dentro de malloc y el proceso moria a mitad de
 * escritura en vez de emitir el reporte.
 *
 * Bloqueando la señal durante el solve, queda pendiente y se entrega al
 * restaurar la mascara -- entre solves, con el heap consistente. El costo es
 * un sobrepaso de a lo sumo un solve (~8 s sobre 1080 s: 0,7%).
 *---------------------------------------------------------------------------*/
namespace {
class AlarmGuard {
public:
    AlarmGuard() {
        sigset_t block;
        sigemptyset(&block);
        sigaddset(&block, SIGALRM);
        sigprocmask(SIG_BLOCK, &block, &prev_);
    }
    ~AlarmGuard() { sigprocmask(SIG_SETMASK, &prev_, nullptr); }
private:
    sigset_t prev_;
};
}  // namespace

FmoResult ImrtFmoSolver::solve(const std::vector<int>& active_angles)
{
    AlarmGuard guard;   // la señal solo se entrega al salir de este scope

    if (!ready_)
        return {std::vector<double>(source_.n_dimlets(), 0.0), 1e30,
                std::vector<double>(ptv_bounds_.size(), 0.0),
                std::vector<double>(oar_bounds_.size(), 0.0)};

    // Instrumentacion de fases: mide cuanto del tiempo por solve es armado de
    // tuplas en C++, cuanto es transferirlas a AMPL, y cuanto es Gurobi. Solo se
    // imprime bajo verbose_; el costo de los relojes es despreciable frente al solve.
    using clk = std::chrono::steady_clock;
    auto t_a = clk::now();

    std::vector<int> active = source_.activeDimletIds(active_angles);

    std::vector<double> active_d(active.begin(), active.end());

    std::vector<ampl::Tuple> ptv_tuples, oar_tuples;
    std::vector<double> ptv_vals, oar_vals;
    std::vector<int> ptv_nnz(ptv_bounds_.size(), 0), oar_nnz(oar_bounds_.size(), 0);

    auto organOf = [](const std::vector<OrganBounds>& bounds, int row) -> int {
        for (size_t i = 0; i < bounds.size(); ++i)
            if (row >= bounds[i].row_off && row < bounds[i].row_off + bounds[i].n_boxets)
                return (int)i;
        return -1;
    };

    for (int j : active) {
        for (const auto& e : source_.ptvDoseFor(j)) {
            ptv_tuples.emplace_back(ampl::Variant((double)e.first), ampl::Variant((double)j));
            ptv_vals.push_back(e.second);
            // Solo alimenta printSolveDims bajo verbose_. organOf es un escaneo
            // lineal por entrada de dosis (~850k por solve), asi que fuera de
            // verbose_ es trabajo que nadie lee.
            if (verbose_) {
                int oi = organOf(ptv_bounds_, e.first);
                if (oi >= 0) ++ptv_nnz[oi];
            }
        }
        for (const auto& e : source_.oarDoseFor(j)) {
            oar_tuples.emplace_back(ampl::Variant((double)e.first), ampl::Variant((double)j));
            oar_vals.push_back(e.second);
            // Solo alimenta printSolveDims bajo verbose_. organOf es un escaneo
            // lineal por entrada de dosis (~850k por solve), asi que fuera de
            // verbose_ es trabajo que nadie lee.
            if (verbose_) {
                int oi = organOf(oar_bounds_, e.first);
                if (oi >= 0) ++oar_nnz[oi];
            }
        }
    }

    auto t_b = clk::now();

    if (verbose_)
        printSolveDims((int)active.size(), ptv_nnz, oar_nnz);

    try {
        // Reset before reassigning: narrowing DIMLETS/PTV_DOSE/OAR_DOSE via `let`
        // while d_ptv/d_oar still hold values for now-removed indices makes AMPL
        // raise "invalid subscripts discarded" as a hard error (not a warning) on
        // that assignment -- which broke every solve() after the first. `reset
        // data` wipes the set and its dependent params atomically with no
        // narrowing step in between, so there's nothing left to discard.
        ampl_->eval("reset data DIMLETS, PTV_DOSE, OAR_DOSE, d_ptv, d_oar;");

        ampl_->getSet("DIMLETS").setValues(ampl::Args(active_d.data()), active_d.size());
        ampl_->getSet("PTV_DOSE").setValues(ptv_tuples.data(), ptv_tuples.size());
        ampl_->getSet("OAR_DOSE").setValues(oar_tuples.data(), oar_tuples.size());

        if (!ptv_tuples.empty())
            ampl_->getParameter("d_ptv").setValues(ptv_tuples.data(), ampl::Args(ptv_vals.data()), ptv_tuples.size());
        if (!oar_tuples.empty())
            ampl_->getParameter("d_oar").setValues(oar_tuples.data(), ampl::Args(oar_vals.data()), oar_tuples.size());

        auto t_c = clk::now();

        ampl_->solve();

        auto t_d = clk::now();

        // AMPL conserva el objetivo de la instancia ANTERIOR cuando el solve no
        // se completa -- licencia caida, modelo infactible, solver interrumpido --
        // y getObjective().value() lo devuelve sin senalar nada. Ese valor viaja
        // como si fuera una evaluacion legitima. En las corridas afectadas se
        // observaron decenas de configuraciones de angulos distintas compartiendo
        // el mismo objetivo hasta el decimal, y busquedas que cortaban por falso
        // minimo local: ningun vecino "mejoraba" porque todos devolvian lo mismo.
        // El mensaje del solver quedaba en el log, pero nada lo miraba.
        //
        // Verificar solve_result convierte ese fallo silencioso en excepcion, que
        // el catch de abajo traduce a 1e30. El vecino queda descartado en vez de
        // contaminar la trayectoria, y la corrida continua.
        const std::string solve_result = ampl_->getValue("solve_result").str();
        if (solve_result != "solved") {
            throw std::runtime_error("solve_result=" + solve_result);
        }

        double f = ampl_->getObjective("fmo_objective").value();

        std::vector<double> x_full(source_.n_dimlets(), 0.0);
        ampl::Variable xvar = ampl_->getVariable("x");
        for (int j : active)
            x_full[j] = xvar.get(ampl::Variant((double)j)).value();

        // Desglose por órgano de lo que el objetivo agregado esconde: suma
        // de u_b^2 (subdosis PTV) y v_b^2 (sobredosis OAR) sobre los boxets
        // de cada órgano. Mismas variables que ya arma fmo.mod, solo que acá
        // se leen y se agregan por órgano en vez de sumarlas todas en una.
        ampl::Variable uvar = ampl_->getVariable("u");
        std::vector<double> ptv_underdose_sq(ptv_bounds_.size(), 0.0);
        for (size_t oi = 0; oi < ptv_bounds_.size(); ++oi) {
            const OrganBounds& b = ptv_bounds_[oi];
            double acc = 0.0;
            for (int row = b.row_off; row < b.row_off + b.n_boxets; ++row) {
                double u = uvar.get(ampl::Variant((double)row)).value();
                acc += u * u;
            }
            ptv_underdose_sq[oi] = acc;
        }

        ampl::Variable vvar = ampl_->getVariable("v");
        std::vector<double> oar_overdose_sq(oar_bounds_.size(), 0.0);
        for (size_t oi = 0; oi < oar_bounds_.size(); ++oi) {
            const OrganBounds& b = oar_bounds_[oi];
            double acc = 0.0;
            for (int row = b.row_off; row < b.row_off + b.n_boxets; ++row) {
                double v = vvar.get(ampl::Variant((double)row)).value();
                acc += v * v;
            }
            oar_overdose_sq[oi] = acc;
        }

        if (verbose_) {
            auto ms = [](clk::time_point a, clk::time_point b) {
                return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
            };
            auto t_e = clk::now();
            std::cout << "    [FMO fases ms] armado=" << ms(t_a, t_b)
                      << " envio=" << ms(t_b, t_c)
                      << " gurobi=" << ms(t_c, t_d)
                      << " lectura=" << ms(t_d, t_e)
                      << " total=" << ms(t_a, t_e) << "\n";
        }

        return {x_full, f, ptv_underdose_sq, oar_overdose_sq};
    } catch (const std::exception& e) {
        std::cerr << "[FMO] ERROR: AMPL/Gurobi solve failed: " << e.what() << "\n";
        return {std::vector<double>(source_.n_dimlets(), 0.0), 1e30,
                std::vector<double>(ptv_bounds_.size(), 0.0),
                std::vector<double>(oar_bounds_.size(), 0.0)};
    }
}

} // namespace imrt
} // namespace emili
