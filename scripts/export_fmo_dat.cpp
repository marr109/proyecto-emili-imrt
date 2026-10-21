// One-off export tool: dumps the exact FMO data (params + sets) that
// ImrtFmoSolver::solve() pushes into AMPL for a given active-angle set, as a
// standalone .dat file. Lets an outside reader run `ampl_gurobi/fmo.mod`
// with a plain `ampl` session and reproduce our objective value, without
// needing our C++/AMPL pipeline. Mirrors imrt_fmo.cpp::solve() exactly
// (same tuple layout, same PTV_B/OAR_B row order) -- not wired into
// CMakeLists.txt, compile ad hoc:
//
//   clang++ -std=c++17 -I. -o /tmp/export_fmo_dat \
//       scripts/export_fmo_dat.cpp imrt/cerr_instance.cpp
//
// Usage:
//   export_fmo_dat <instance_dir> <output.dat> <angle_deg>...

#include "imrt/cerr_instance.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using emili::imrt::CerrFmoSource;
using emili::imrt::FmoOrganRef;

namespace {

bool allEqual(const std::vector<double>& v, double eps = 1e-9)
{
    for (double x : v)
        if (std::fabs(x - v[0]) > eps) return false;
    return true;
}

// Writes `param name := v1 v2 ...;` (1 value per PTV_B/OAR_B row, in the
// same concatenated organ order ptvOrgans()/oarOrgans() build), or the
// compact `param name default v;` form when every row shares one value --
// both are semantically identical for AMPL, the second is just far shorter
// when true (as it is for every clinical value in this placeholder config).
// allow_default must be false for params that already declare a `default`
// in fmo.mod itself (dmax_ptv) -- AMPL rejects a second default from data
// ("default overrides that in model"), so those need the explicit list.
void writeParamVector(std::ofstream& out, const char* name, const std::vector<double>& v,
                       bool allow_default = true)
{
    if (v.empty()) return;
    if (allow_default && allEqual(v)) {
        out << "param " << name << " default " << v[0] << ";\n";
        return;
    }
    out << "param " << name << " :=\n";
    for (size_t i = 0; i < v.size(); ++i)
        out << "  " << i << " " << v[i] << "\n";
    out << ";\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "usage: export_fmo_dat <instance_dir> <output.dat> <angle_deg>...\n";
        return 1;
    }
    std::string instance_dir = argv[1];
    std::string output_path  = argv[2];
    std::vector<int> active_angles;
    for (int i = 3; i < argc; ++i) active_angles.push_back(std::atoi(argv[i]));

    CerrFmoSource src(instance_dir);

    // Same accumulation order as ImrtFmoSolver::precompute() (imrt_fmo.cpp).
    int n_ptv = 0, n_oar = 0;
    std::vector<double> dmin, dmax_ptv, dmax;
    for (const FmoOrganRef& o : src.ptvOrgans()) {
        n_ptv += o.n_boxets;
        for (int b = 0; b < o.n_boxets; ++b) {
            dmin.push_back(o.dmin);
            dmax_ptv.push_back(o.dmax_ptv);
        }
    }
    for (const FmoOrganRef& o : src.oarOrgans()) {
        n_oar += o.n_boxets;
        for (int b = 0; b < o.n_boxets; ++b)
            dmax.push_back(o.dmax);
    }

    std::vector<int> active = src.activeDimletIds(active_angles);

    std::ofstream out(output_path);
    out << std::setprecision(12);

    out << "# Auto-generado por scripts/export_fmo_dat.cpp -- refleja exactamente\n"
        << "# los datos que ImrtFmoSolver::solve() envia a AMPL para este set de\n"
        << "# angulos activos. Usar junto con fmo.mod:\n"
        << "#   ampl: model fmo.mod;\n"
        << "#   ampl: data " << output_path << ";\n"
        << "#   ampl: option solver gurobi; solve; display fmo_objective;\n"
        << "# Angulos activos (grados): ";
    for (int a : active_angles) out << a << " ";
    out << "\n\n";

    out << "param n_ptv := " << n_ptv << ";\n";
    out << "param n_oar := " << n_oar << ";\n";
    out << "param max_intensity := " << src.max_intensity() << ";\n";
    out << "param w_under := " << src.w_under() << ";\n";
    out << "param w_over := " << src.w_over() << ";\n";
    out << "param w_ptv_over := " << src.w_ptv_over() << ";\n\n";

    writeParamVector(out, "dmin", dmin);
    writeParamVector(out, "dmax_ptv", dmax_ptv, /*allow_default=*/false);
    writeParamVector(out, "dmax", dmax);
    out << "\n";

    out << "set DIMLETS :=\n  ";
    for (size_t i = 0; i < active.size(); ++i) out << active[i] << " ";
    out << "\n;\n\n";

    // Same loop as ImrtFmoSolver::solve() (imrt_fmo.cpp): for each active
    // dimlet, pull its sparse (boxet_row, dose) hits against PTV and OAR.
    std::vector<std::pair<int,int>>   ptv_pairs; // (boxet, dimlet)
    std::vector<double>               ptv_vals;
    std::vector<std::pair<int,int>>   oar_pairs;
    std::vector<double>               oar_vals;

    for (int j : active) {
        for (const auto& e : src.ptvDoseFor(j)) {
            ptv_pairs.emplace_back(e.first, j);
            ptv_vals.push_back(e.second);
        }
        for (const auto& e : src.oarDoseFor(j)) {
            oar_pairs.emplace_back(e.first, j);
            oar_vals.push_back(e.second);
        }
    }

    out << "set PTV_DOSE :=\n";
    for (const auto& p : ptv_pairs) out << "  (" << p.first << "," << p.second << ")\n";
    out << ";\n\n";

    out << "param d_ptv :=\n";
    for (size_t i = 0; i < ptv_pairs.size(); ++i)
        out << "  [" << ptv_pairs[i].first << "," << ptv_pairs[i].second << "] " << ptv_vals[i] << "\n";
    out << ";\n\n";

    out << "set OAR_DOSE :=\n";
    for (const auto& p : oar_pairs) out << "  (" << p.first << "," << p.second << ")\n";
    out << ";\n\n";

    out << "param d_oar :=\n";
    for (size_t i = 0; i < oar_pairs.size(); ++i)
        out << "  [" << oar_pairs[i].first << "," << oar_pairs[i].second << "] " << oar_vals[i] << "\n";
    out << ";\n";

    std::cerr << "wrote " << output_path
              << " (n_ptv=" << n_ptv << " n_oar=" << n_oar
              << " active_dimlets=" << active.size()
              << " ptv_entries=" << ptv_pairs.size()
              << " oar_entries=" << oar_pairs.size() << ")\n";
    return 0;
}
