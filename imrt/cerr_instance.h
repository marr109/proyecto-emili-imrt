#ifndef CERR_INSTANCE_H
#define CERR_INSTANCE_H

#include "imrt_fmo_source.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>

namespace emili {
namespace imrt {

/*---------------------------------------------------------------------------*
 * CerrFmoSource
 *
 * Native C++ reader for the raw CERR export layout in instances/CERR_Prostate.
 * Angle beamlet counts are read directly from beamletIndex.txt
 * instead of assumed uniform, since this dataset does not have a fixed
 * beamlet count per angle.
 *
 * Layout read:
 *   <ORGAN>.txt              global voxel ids; line order = local boxet index
 *   <ORGAN>_<angle_idx>.txt  rows: global_voxel_id local_beamlet_id(1-based) dose_rate
 *   beamletIndex.txt         rows: angle_idx global_start global_end (1-based, inclusive)
 *
 * The per-angle dose files are ~11GB total across 360 angles x 4 organs, so
 * they are parsed lazily, one angle at a time, on first request -- BAO and
 * the classic imrt local search only touch a subset of angles per move, and
 * revisit angles across iterations, so caching the parsed rows in memory
 * (unbounded, no eviction) avoids redundant disk parsing without needing a
 * size cap.
 *---------------------------------------------------------------------------*/
class CerrFmoSource : public IFmoDataSource {
public:
    explicit CerrFmoSource(const std::string& dir);

    int    n_dimlets()     const override { return n_dimlets_total_; }
    double max_intensity() const override { return max_intensity_; }
    double w_under()       const override { return w_under_; }
    double w_over()        const override { return w_over_; }
    double w_ptv_over()    const override { return w_ptv_over_; }

    const std::vector<FmoOrganRef>& ptvOrgans() const override { return ptv_organs_; }
    const std::vector<FmoOrganRef>& oarOrgans() const override { return oar_organs_; }

    std::vector<int> activeDimletIds(const std::vector<int>& active_angles) const override;

    const std::vector<std::pair<int,double>>& ptvDoseFor(int global_dimlet_id) const override;
    const std::vector<std::pair<int,double>>& oarDoseFor(int global_dimlet_id) const override;

    int nAnglesTotal() const { return n_angles_total_; }

    // True if `dir` looks like a CERR-format export (beamletIndex.txt present,
    // instance_config.txt absent) -- used by ImrtBuilder to auto-detect the
    // instance format before constructing this source.
    static bool looksLikeCerrDir(const std::string& dir);

private:
    struct OrganMeta {
        std::string name;
        int row_off;   // offset of this organ's boxets in the PTV/OAR block
        int n_boxets;
        std::unordered_map<int,int> voxel_to_local;
    };

    std::string dir_;
    int n_angles_total_;
    int n_dimlets_total_;
    std::vector<std::pair<int,int>> beamlet_range_; // angle_idx -> (start,end), 1-based inclusive

    std::vector<OrganMeta> ptv_meta_;
    std::vector<OrganMeta> oar_meta_;
    std::vector<FmoOrganRef> ptv_organs_;
    std::vector<FmoOrganRef> oar_organs_;

    double max_intensity_, w_under_, w_over_, w_ptv_over_;

    mutable std::vector<bool> angle_loaded_;
    mutable std::vector<std::vector<std::pair<int,double>>> ptv_dose_;
    mutable std::vector<std::vector<std::pair<int,double>>> oar_dose_;

    void loadBeamletIndex();
    // skip: voxel ids to leave out (OAR voxels that also belong to the PTV).
    void loadVoxelList(const std::string& organ_name, OrganMeta& out,
                       const std::unordered_set<int>* skip = nullptr) const;
    void ensureAngleLoaded(int angle_idx) const;
    int  angleForDimlet(int global_dimlet_id) const;
};

} // namespace imrt
} // namespace emili

#endif // CERR_INSTANCE_H
