#ifndef IMRT_FMO_SOURCE_H
#define IMRT_FMO_SOURCE_H

#include <string>
#include <utility>
#include <vector>

namespace emili {
namespace imrt {

/*---------------------------------------------------------------------------*
 * IFmoDataSource
 *
 * Everything ImrtFmoSolver/ImrtProblem needs to build the FMO QP or evaluate
 * a beamlet intensity vector, decoupled from the concrete on-disk instance
 * format. CerrFmoSource (instances/CERR_Prostate) has a variable beamlet
 * count per angle, so this interface avoids assuming a fixed per-angle
 * stride the way the old ImrtInstance::globalDimletIndex() did.
 *---------------------------------------------------------------------------*/
struct FmoOrganRef {
    std::string name;
    bool   is_ptv;
    int    n_boxets;
    double dmin;       // PTV prescription (Gy); unused for OAR
    double dmax;       // OAR tolerance (Gy); unused for PTV
    double dmax_ptv;   // PTV overdose ceiling (e.g. 1.07*dmin); unused for OAR
};

class IFmoDataSource {
public:
    virtual ~IFmoDataSource() {}

    virtual int    n_dimlets()     const = 0;
    virtual double max_intensity() const = 0;
    virtual double w_under()       const = 0;
    virtual double w_over()        const = 0;
    virtual double w_ptv_over()    const = 0;

    // In boxet-block order: row r of ptvDoseFor()/oarDoseFor() falls within
    // the r-th organ's [0, n_boxets) slice of these lists, concatenated.
    virtual const std::vector<FmoOrganRef>& ptvOrgans() const = 0;
    virtual const std::vector<FmoOrganRef>& oarOrgans() const = 0;

    // Global dimlet ids that are active for this set of active angle indices.
    virtual std::vector<int> activeDimletIds(const std::vector<int>& active_angles) const = 0;

    // Sparse dose entries touching a given global dimlet: (boxet row within the
    // concatenated PTV/OAR block, dose_rate). Empty if the dimlet touches nothing.
    virtual const std::vector<std::pair<int,double>>& ptvDoseFor(int global_dimlet_id) const = 0;
    virtual const std::vector<std::pair<int,double>>& oarDoseFor(int global_dimlet_id) const = 0;

    // Computes per-boxet-row PTV/OAR doses for a full intensity vector x
    // (length n_dimlets()), in the same concatenated row order as
    // ptvOrgans()/oarOrgans(). Skips zero-intensity dimlets so a search that
    // never touches most of a lazily-loaded source (CerrFmoSource) doesn't
    // force every angle's dose file to be parsed just to evaluate a solution.
    void computeDoses(const std::vector<double>& x,
                       std::vector<double>& ptv_dose,
                       std::vector<double>& oar_dose) const;
};

} // namespace imrt
} // namespace emili

#endif // IMRT_FMO_SOURCE_H
