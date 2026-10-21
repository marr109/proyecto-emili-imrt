#include "cerr_instance.h"
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <tuple>

namespace emili {
namespace imrt {

namespace {
const char* kPtvNames[] = { "PTVHD", "PTVLD" };
const char* kOarNames[] = { "BLADDER", "RECTUM" };

std::string withTrailingSlash(const std::string& dir)
{
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') return dir + "/";
    return dir;
}
} // namespace

bool CerrFmoSource::looksLikeCerrDir(const std::string& dir)
{
    std::string d = withTrailingSlash(dir);
    std::ifstream beamlet_index(d + "beamletIndex.txt");
    std::ifstream config(d + "instance_config.txt");
    return beamlet_index.good() && !config.good();
}

CerrFmoSource::CerrFmoSource(const std::string& dir)
    : dir_(withTrailingSlash(dir))
    , n_angles_total_(0), n_dimlets_total_(0)
{
    loadBeamletIndex();

    angle_loaded_.assign(n_angles_total_, false);
    ptv_dose_.assign(n_dimlets_total_, {});
    oar_dose_.assign(n_dimlets_total_, {});

    // Clinical Dmin/Dmax placeholders -- authorized for pipeline validation
    // only (no instance_config.txt / CERR patient plan documents the real
    // prescription for this dataset).
    // EMILI_DMIN overrides the PTV prescription (Gy) for both PTVs; dmax_ptv
    // follows as 1.07 * dmin. Objectives are not comparable across values.
    double dmin_rx = 65.0;
    if (const char* env_d = std::getenv("EMILI_DMIN")) dmin_rx = std::atof(env_d);
    const double dmin_ptvhd = dmin_rx, dmin_ptvld = dmin_rx;
    const double dmax_bladder = 50.0, dmax_rectum = 50.0;
    // EMILI_MAX_INTENSITY overrides the per-beamlet intensity cap. The default (15000) is the
    // original placeholder and never binds; a clinically meaningful cap is set per experiment.
    max_intensity_ = 15000.0;
    if (const char* env_mi = std::getenv("EMILI_MAX_INTENSITY")) {
        max_intensity_ = std::atof(env_mi);
    }
    w_under_       = 1.0;
    w_over_        = 0.5;
    // La penalizacion de sobredosis en PTV esta DESACTIVADA por defecto. No es
    // un olvido: el modelo optimiza cobertura del PTV y proteccion de los OAR,
    // no homogeneidad, y con w_ptv_over_ = 0 el termino w de fmo.mod queda sin
    // costo y la fila ptv_ceiling no ata. La consecuencia es que la dosis en el
    // PTV puede crecer sin penalizacion --Dmax observado ~5x la prescripcion--
    // y los indices de conformidad y homogeneidad no describen calidad de plan.
    //
    // EMILI_W_PTV_OVER activa el termino sin recompilar, para poder medir su
    // efecto. ATENCION: al activarlo el objetivo deja de ser comparable con las
    // corridas previas. No es el mismo problema mejor resuelto, es otro problema.
    w_ptv_over_    = 0.0;
    if (const char* env_w = std::getenv("EMILI_W_PTV_OVER")) {
        w_ptv_over_ = std::atof(env_w);
    }

    int off = 0;
    for (const char* name : kPtvNames) {
        OrganMeta m;
        m.name = name;
        m.row_off = off;
        loadVoxelList(name, m);
        off += m.n_boxets;

        double dmin = (m.name == "PTVHD") ? dmin_ptvhd : dmin_ptvld;
        FmoOrganRef ref;
        ref.name     = name;
        ref.is_ptv   = true;
        ref.n_boxets = m.n_boxets;
        ref.dmin     = dmin;
        ref.dmax     = 0.0;
        ref.dmax_ptv = 1.07 * dmin;
        ptv_organs_.push_back(ref);
        ptv_meta_.push_back(std::move(m));
    }

    // EMILI_EXCLUDE_OVERLAP=1 -> OAR minus PTV: voxels that belong to both a PTV
    // and an OAR count only as PTV. Without it the same voxel must satisfy
    // dose >= dmin (PTV) and dose <= dmax (OAR), which is contradictory.
    // Changes the problem: objectives are not comparable with runs without it.
    std::unordered_set<int> ptv_voxel_ids;
    const char* env_ov = std::getenv("EMILI_EXCLUDE_OVERLAP");
    const bool exclude_overlap = env_ov && std::string(env_ov) == "1";
    if (exclude_overlap)
        for (const auto& pm : ptv_meta_)
            for (const auto& kv : pm.voxel_to_local) ptv_voxel_ids.insert(kv.first);

    off = 0;
    for (const char* name : kOarNames) {
        OrganMeta m;
        m.name = name;
        m.row_off = off;
        loadVoxelList(name, m, exclude_overlap ? &ptv_voxel_ids : nullptr);
        off += m.n_boxets;

        double dmax = (m.name == "BLADDER") ? dmax_bladder : dmax_rectum;
        FmoOrganRef ref;
        ref.name     = name;
        ref.is_ptv   = false;
        ref.n_boxets = m.n_boxets;
        ref.dmin     = 0.0;
        ref.dmax     = dmax;
        ref.dmax_ptv = 0.0;
        oar_organs_.push_back(ref);
        oar_meta_.push_back(std::move(m));
    }
}

void CerrFmoSource::loadBeamletIndex()
{
    std::ifstream f(dir_ + "beamletIndex.txt");
    std::vector<std::tuple<int,int,int>> rows;
    std::string line;
    int max_idx = -1, max_end = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int idx, s, e;
        if (!(iss >> idx >> s >> e)) continue;
        rows.emplace_back(idx, s, e);
        if (idx > max_idx) max_idx = idx;
        if (e > max_end)   max_end = e;
    }

    n_angles_total_  = max_idx + 1;
    n_dimlets_total_ = max_end;
    beamlet_range_.assign(n_angles_total_, {0, 0});
    for (const auto& row : rows)
        beamlet_range_[std::get<0>(row)] = {std::get<1>(row), std::get<2>(row)};
}

void CerrFmoSource::loadVoxelList(const std::string& organ_name, OrganMeta& out,
                                  const std::unordered_set<int>* skip) const
{
    std::ifstream f(dir_ + organ_name + ".txt");
    std::string line;
    int local = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int voxel_id;
        if (!(iss >> voxel_id)) continue;
        if (skip && skip->count(voxel_id)) continue;
        out.voxel_to_local[voxel_id] = local++;
    }
    out.n_boxets = local;
}

void CerrFmoSource::ensureAngleLoaded(int angle_idx) const
{
    if (angle_idx < 0 || angle_idx >= n_angles_total_) return;
    if (angle_loaded_[angle_idx]) return;

    const int g_start = beamlet_range_[angle_idx].first;

    auto loadOrgan = [&](const OrganMeta& m, std::vector<std::vector<std::pair<int,double>>>& dose) {
        std::ifstream f(dir_ + m.name + "_" + std::to_string(angle_idx) + ".txt");
        if (!f.is_open()) return;

        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            int voxel_id, local_beamlet;
            double dose_rate;
            if (!(iss >> voxel_id >> local_beamlet >> dose_rate)) continue;

            auto it = m.voxel_to_local.find(voxel_id);
            if (it == m.voxel_to_local.end()) continue;

            int global_dimlet = (g_start - 1) + (local_beamlet - 1);
            if (global_dimlet < 0 || global_dimlet >= n_dimlets_total_) continue;

            dose[global_dimlet].push_back({m.row_off + it->second, dose_rate});
        }
    };

    for (const auto& m : ptv_meta_) loadOrgan(m, ptv_dose_);
    for (const auto& m : oar_meta_) loadOrgan(m, oar_dose_);

    angle_loaded_[angle_idx] = true;
}

int CerrFmoSource::angleForDimlet(int global_dimlet_id) const
{
    // beamlet_range_ is contiguous and sorted by angle_idx (0..n_angles_total_-1),
    // so binary search on the start bound locates the owning angle.
    int g = global_dimlet_id + 1;
    int lo = 0, hi = n_angles_total_ - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (beamlet_range_[mid].first <= g) lo = mid; else hi = mid - 1;
    }
    return lo;
}

std::vector<int> CerrFmoSource::activeDimletIds(const std::vector<int>& active_angles) const
{
    std::vector<int> ids;
    for (int a : active_angles) {
        ensureAngleLoaded(a);
        if (a < 0 || a >= n_angles_total_) continue;
        const auto& r = beamlet_range_[a];
        for (int g = r.first - 1; g <= r.second - 1; ++g) ids.push_back(g);
    }
    return ids;
}

const std::vector<std::pair<int,double>>& CerrFmoSource::ptvDoseFor(int global_dimlet_id) const
{
    ensureAngleLoaded(angleForDimlet(global_dimlet_id));
    return ptv_dose_[global_dimlet_id];
}

const std::vector<std::pair<int,double>>& CerrFmoSource::oarDoseFor(int global_dimlet_id) const
{
    ensureAngleLoaded(angleForDimlet(global_dimlet_id));
    return oar_dose_[global_dimlet_id];
}

} // namespace imrt
} // namespace emili
