#include "imrt_fmo_source.h"
#include <algorithm>

namespace emili {
namespace imrt {

void IFmoDataSource::computeDoses(const std::vector<double>& x,
                                   std::vector<double>& ptv_dose,
                                   std::vector<double>& oar_dose) const
{
    int n_ptv = 0, n_oar = 0;
    for (const FmoOrganRef& o : ptvOrgans()) n_ptv += o.n_boxets;
    for (const FmoOrganRef& o : oarOrgans()) n_oar += o.n_boxets;

    ptv_dose.assign(n_ptv, 0.0);
    oar_dose.assign(n_oar, 0.0);

    int n = n_dimlets();
    int lim = std::min(n, (int)x.size());
    for (int d = 0; d < lim; ++d) {
        if (x[d] == 0.0) continue;
        for (const auto& e : ptvDoseFor(d))
            ptv_dose[e.first] += e.second * x[d];
        for (const auto& e : oarDoseFor(d))
            oar_dose[e.first] += e.second * x[d];
    }
}

} // namespace imrt
} // namespace emili
