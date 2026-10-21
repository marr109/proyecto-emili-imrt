#include "imrt.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace emili {
namespace imrt {

/*---------------------------------------------------------------------------*
 * ImrtProblem
 *---------------------------------------------------------------------------*/

// Solo se penalizan las violaciones (cuadrático, de un solo lado): un boxet
// PTV en o por encima de Dmin, o un boxet OAR en o por debajo de Dmax, aporta 0.
// Ejemplo: Dmin=65 Gy, un boxet PTV con dosis=60 Gy -> under=5 -> +w_under*25 a f.
double ImrtProblem::calcObjectiveFunctionValue(emili::Solution &solution) {
  ImrtSolution &sol = static_cast<ImrtSolution &>(solution);
  const std::vector<double> &x = sol.getIntensities();

  std::vector<double> ptv_dose, oar_dose;
  source_->computeDoses(x, ptv_dose, oar_dose);

  double f = 0.0;
  double w_under = source_->w_under();
  double w_over  = source_->w_over();

  int row = 0;
  for (const FmoOrganRef &org : source_->ptvOrgans()) {
    for (int b = 0; b < org.n_boxets; ++b, ++row) {
      double under = org.dmin - ptv_dose[row];
      if (under > 0.0)
        f += w_under * under * under;
    }
  }

  row = 0;
  for (const FmoOrganRef &org : source_->oarOrgans()) {
    for (int b = 0; b < org.n_boxets; ++b, ++row) {
      double over = oar_dose[row] - org.dmax;
      if (over > 0.0)
        f += w_over * over * over;
    }
  }

  return f;
}

double ImrtProblem::evaluateSolution(emili::Solution &solution) {
  double val = calcObjectiveFunctionValue(solution);
  solution.setSolutionValue(val);
  if (verbose_fmo_ && val < last_best_fmo_) {
    double delta = last_best_fmo_ - val; // positive = improvement
    printFmoStats(static_cast<ImrtSolution &>(solution), val, delta);
    last_best_fmo_ = val;
  }
  return val;
}

void ImrtProblem::setActiveAngles(int k) {
  active_angle_idxs_.clear();
  k = std::min(k, nAngles());
  for (int i = 0; i < k; ++i)
    active_angle_idxs_.push_back(i);
}

bool ImrtProblem::isAngleActive(int angle_idx) const {
  if (active_angle_idxs_.empty())
    return true;
  for (int ai : active_angle_idxs_)
    if (ai == angle_idx)
      return true;
  return false;
}

void ImrtProblem::printFmoStats(const ImrtSolution &sol, double val,
                                double delta) {
  const auto &x = sol.getIntensities();
  std::vector<double> ptv_dose, oar_dose;
  source_->computeDoses(x, ptv_dose, oar_dose);

  int n_angles = nAngles();
  int n_active = active_angle_idxs_.empty() ? n_angles
                                             : (int)active_angle_idxs_.size();

  std::cout << "\n------------------------------------------------------------\n";
  ++improve_count_;
  std::cout << " FMO mejora #" << improve_count_
            << " | f = " << std::fixed << std::setprecision(2) << val;
  if (improve_count_ > 1)
    std::cout << " | Delta = -" << std::setprecision(2) << delta;
  std::cout << "\n------------------------------------------------------------\n";

  // Per-angle stats
  std::cout << " Angulos gantry candidatos (" << n_active << "/"
            << n_angles << " activos):\n";
  for (int a = 0; a < n_angles; ++a) {
    if (isAngleActive(a)) {
      std::vector<int> ids = source_->activeDimletIds({a});
      double sum_x = 0.0;
      int nonzero = 0;
      for (int d : ids) {
        sum_x += x[d];
        if (x[d] > 0.0)
          ++nonzero;
      }
      int per = (int)ids.size();
      std::cout << "   [ACT] " << std::setw(4) << angle_degrees_[a]
                << " deg | no-cero: " << std::setw(3) << nonzero << "/" << per
                << " | media: " << std::setprecision(3)
                << (per > 0 ? sum_x / per : 0.0)
                << " | suma: " << std::setprecision(2) << sum_x << "\n";
    } else {
      std::cout << "   [---] " << std::setw(4) << angle_degrees_[a]
                << " deg | (inactivo)\n";
    }
  }

  // Per-organ stats
  std::cout << " Organos:\n";
  int row = 0;
  for (const FmoOrganRef &org : source_->ptvOrgans()) {
    double mean_d = 0.0, penalty = 0.0;
    for (int b = 0; b < org.n_boxets; ++b) {
      double d = ptv_dose[row + b];
      mean_d += d;
      double u = org.dmin - d;
      if (u > 0)
        penalty += source_->w_under() * u * u;
    }
    if (org.n_boxets > 0)
      mean_d /= (double)org.n_boxets;
    row += org.n_boxets;

    std::cout << "   " << std::left << std::setw(10) << org.name
              << " [tumor] | dosis media: " << std::setprecision(2)
              << mean_d << " Gy | objetivo >= " << org.dmin
              << " | penaliz: " << std::setprecision(2) << penalty << "\n";
  }

  row = 0;
  for (const FmoOrganRef &org : source_->oarOrgans()) {
    double mean_d = 0.0, penalty = 0.0;
    for (int b = 0; b < org.n_boxets; ++b) {
      double d = oar_dose[row + b];
      mean_d += d;
      double ov = d - org.dmax;
      if (ov > 0)
        penalty += source_->w_over() * ov * ov;
    }
    if (org.n_boxets > 0)
      mean_d /= (double)org.n_boxets;
    row += org.n_boxets;

    std::cout << "   " << std::left << std::setw(10) << org.name
              << " [OAR]   | dosis media: " << std::setprecision(2)
              << mean_d << " Gy | limite   <= " << org.dmax
              << " | penaliz: " << std::setprecision(2) << penalty << "\n";
  }
  std::cout << "------------------------------------------------------------\n";
  std::cout << std::right; // restore default alignment
}

/*---------------------------------------------------------------------------*
 * ImrtSolution
 *---------------------------------------------------------------------------*/

const void *ImrtSolution::getRawData() const {
  return static_cast<const void *>(&intensities_);
}

void ImrtSolution::setRawData(const void *data) {
  intensities_ = *static_cast<const std::vector<double> *>(data);
}

std::string ImrtSolution::getSolutionRepresentation() {
  std::ostringstream oss;
  oss << "[";
  for (int i = 0; i < static_cast<int>(intensities_.size()); ++i) {
    if (i > 0)
      oss << ", ";
    oss << intensities_[i];
  }
  oss << "]";
  return oss.str();
}

emili::Solution *ImrtSolution::clone() {
  return new ImrtSolution(solution_value, intensities_);
}

bool ImrtSolution::isFeasible() {
  for (double v : intensities_)
    if (v < 0.0)
      return false;
  return true;
}

/*---------------------------------------------------------------------------*
 * Initial solutions
 *---------------------------------------------------------------------------*/

emili::Solution *ImrtInitialSolution::generateSolution() {
  emili::Solution *s = generate();
  problem_.evaluateSolution(*s);
  return s;
}

emili::Solution *ImrtInitialSolution::generateEmptySolution() {
  return new ImrtSolution(problem_.getSource().n_dimlets());
}

emili::Solution *ZeroInitialSolution::generate() {
  return new ImrtSolution(problem_.getSource().n_dimlets());
}

emili::Solution *UniformInitialSolution::generate() {
  int n = problem_.getSource().n_dimlets();
  const auto &active = problem_.getActiveAngles();
  std::vector<double> x(n, 0.0);
  if (active.empty()) {
    std::fill(x.begin(), x.end(), intensity_);
  } else {
    for (int d : problem_.getSource().activeDimletIds(active))
      x[d] = intensity_;
  }
  return new ImrtSolution(0.0, x);
}

emili::Solution *RandomInitialSolution::generate() {
  int n = problem_.getSource().n_dimlets();
  const auto &active = problem_.getActiveAngles();
  std::vector<double> x(n, 0.0);
  if (active.empty()) {
    for (int d = 0; d < n; ++d)
      x[d] = emili::generateRealRandomNumber() * max_intensity_;
  } else {
    for (int d : problem_.getSource().activeDimletIds(active))
      x[d] = emili::generateRealRandomNumber() * max_intensity_;
  }
  return new ImrtSolution(0.0, x);
}

/*---------------------------------------------------------------------------*
 * ImrtNeighborhood
 *---------------------------------------------------------------------------*/

ImrtNeighborhood::ImrtNeighborhood(ImrtProblem &p)
    : problem_(p), n_dimlets_(p.getSource().n_dimlets()),
      max_intensity_(p.getSource().max_intensity()) {
  if (!p.getActiveAngles().empty())
    buildActiveDimlets(p.getActiveAngles());
}

void ImrtNeighborhood::buildActiveDimlets(const std::vector<int> &angle_idxs) {
  active_dimlets_ = problem_.getSource().activeDimletIds(angle_idxs);
}

/*---------------------------------------------------------------------------*
 * SingleBeamletShift
 *---------------------------------------------------------------------------*/

emili::Neighborhood::NeighborhoodIterator
SingleBeamletShift::begin(emili::Solution *base) {
  cur_dimlet_     = 0;
  cur_active_idx_ = 0;
  cur_dir_        = 1;
  first_          = true;
  return emili::Neighborhood::NeighborhoodIterator(this, base);
}

void SingleBeamletShift::reset() {
  cur_dimlet_     = 0;
  cur_active_idx_ = 0;
  cur_dir_        = 1;
  first_          = true;
}

emili::Solution *SingleBeamletShift::computeStep(emili::Solution *step) {
  if (!active_dimlets_.empty()) {
    if (!first_) {
      if (cur_dir_ == 1) {
        cur_dir_ = -1;
      } else {
        cur_dir_ = 1;
        ++cur_active_idx_;
      }
    }
    first_ = false;
    if (cur_active_idx_ >= (int)active_dimlets_.size())
      return nullptr;
    cur_dimlet_ = active_dimlets_[cur_active_idx_];
  } else {
    if (!first_) {
      if (cur_dir_ == 1) {
        cur_dir_ = -1;
      } else {
        cur_dir_ = 1;
        ++cur_dimlet_;
      }
    }
    first_ = false;
    if (cur_dimlet_ >= n_dimlets_)
      return nullptr;
  }

  ImrtSolution *s = static_cast<ImrtSolution *>(step);
  std::vector<double> x = s->getIntensities();

  last_dimlet_    = cur_dimlet_;
  last_old_value_ = x[cur_dimlet_];

  double newval = x[cur_dimlet_] + cur_dir_ * delta_;
  newval = std::max(0.0, std::min(newval, max_intensity_));
  x[cur_dimlet_] = newval;

  s->setIntensities(x);
  problem_.evaluateSolution(*s);
  return s;
}

void SingleBeamletShift::reverseLastMove(emili::Solution *step) {
  ImrtSolution *s = static_cast<ImrtSolution *>(step);
  std::vector<double> x = s->getIntensities();
  x[last_dimlet_] = last_old_value_;
  s->setIntensities(x);
}

emili::Solution *SingleBeamletShift::random(emili::Solution *s) {
  ImrtSolution *sol = static_cast<ImrtSolution *>(s);
  std::vector<double> x = sol->getIntensities();

  int b;
  if (!active_dimlets_.empty()) {
    b = active_dimlets_[emili::generateRandomNumber() % active_dimlets_.size()];
  } else {
    b = emili::generateRandomNumber() % n_dimlets_;
  }
  int dir = (emili::generateRandomNumber() % 2 == 0) ? 1 : -1;
  x[b] = std::max(0.0, std::min(x[b] + dir * delta_, max_intensity_));

  ImrtSolution *nb = new ImrtSolution(0.0, x);
  problem_.evaluateSolution(*nb);
  return nb;
}

/*---------------------------------------------------------------------------*
 * BeamletSwap
 *---------------------------------------------------------------------------*/

emili::Neighborhood::NeighborhoodIterator
BeamletSwap::begin(emili::Solution *base) {
  i_     = 0; j_     = 1;
  i_act_ = 0; j_act_ = 1;
  first_ = true;
  return emili::Neighborhood::NeighborhoodIterator(this, base);
}

void BeamletSwap::reset() {
  i_     = 0; j_     = 1;
  i_act_ = 0; j_act_ = 1;
  first_ = true;
}

emili::Solution *BeamletSwap::computeStep(emili::Solution *step) {
  if (!active_dimlets_.empty()) {
    int na = (int)active_dimlets_.size();
    if (!first_) {
      ++j_act_;
      if (j_act_ >= na) {
        ++i_act_;
        j_act_ = i_act_ + 1;
      }
    }
    first_ = false;
    if (i_act_ >= na - 1)
      return nullptr;
    last_i_ = active_dimlets_[i_act_];
    last_j_ = active_dimlets_[j_act_];
  } else {
    if (!first_) {
      ++j_;
      if (j_ >= n_dimlets_) {
        ++i_;
        j_ = i_ + 1;
      }
    }
    first_ = false;
    if (i_ >= n_dimlets_ - 1)
      return nullptr;
    last_i_ = i_;
    last_j_ = j_;
  }

  ImrtSolution *s = static_cast<ImrtSolution *>(step);
  std::vector<double> x = s->getIntensities();
  std::swap(x[last_i_], x[last_j_]);

  s->setIntensities(x);
  problem_.evaluateSolution(*s);
  return s;
}

void BeamletSwap::reverseLastMove(emili::Solution *step) {
  ImrtSolution *s = static_cast<ImrtSolution *>(step);
  std::vector<double> x = s->getIntensities();
  std::swap(x[last_i_], x[last_j_]);
  s->setIntensities(x);
}

emili::Solution *BeamletSwap::random(emili::Solution *s) {
  ImrtSolution *sol = static_cast<ImrtSolution *>(s);
  std::vector<double> x = sol->getIntensities();

  int a, b;
  if (!active_dimlets_.empty()) {
    int na = (int)active_dimlets_.size();
    a = active_dimlets_[emili::generateRandomNumber() % na];
    b = active_dimlets_[emili::generateRandomNumber() % na];
    while (b == a)
      b = active_dimlets_[emili::generateRandomNumber() % na];
  } else {
    a = emili::generateRandomNumber() % n_dimlets_;
    b = emili::generateRandomNumber() % n_dimlets_;
    while (b == a)
      b = emili::generateRandomNumber() % n_dimlets_;
  }
  std::swap(x[a], x[b]);

  ImrtSolution *nb = new ImrtSolution(0.0, x);
  problem_.evaluateSolution(*nb);
  return nb;
}

/*---------------------------------------------------------------------------*
 * Perturbation
 *---------------------------------------------------------------------------*/

emili::Solution *RandomBeamletPerturbation::perturb(emili::Solution *solution) {
  ImrtSolution *s = static_cast<ImrtSolution *>(solution);
  std::vector<double> x = s->getIntensities();
  const auto &active = problem_.getActiveAngles();
  int n = static_cast<int>(x.size());

  std::vector<int> active_dimlets;
  if (!active.empty())
    active_dimlets = problem_.getSource().activeDimletIds(active);

  for (int i = 0; i < k_; ++i) {
    int d;
    if (active_dimlets.empty()) {
      d = emili::generateRandomNumber() % n;
    } else {
      d = active_dimlets[emili::generateRandomNumber() % active_dimlets.size()];
    }
    x[d] = emili::generateRealRandomNumber() * max_intensity_;
  }

  ImrtSolution *p = new ImrtSolution(0.0, x);
  problem_.evaluateSolution(*p);
  return p;
}

/*---------------------------------------------------------------------------*
 * Acceptance
 *---------------------------------------------------------------------------*/

emili::Solution *ImrtImproveAccept::accept(emili::Solution *current,
                                           emili::Solution *candidate) {
  return (*candidate < *current) ? candidate : current;
}

/*---------------------------------------------------------------------------*
 * Termination
 *---------------------------------------------------------------------------*/

bool ImrtMaxIterations::terminate(emili::Solution *c, emili::Solution *n) {
  return ++cur_ >= max_;
}

bool ImrtFeasibleTermination::terminate(emili::Solution *c,
                                        emili::Solution *n) {
  return (n != nullptr && n->getSolutionValue() <= tolerance_);
}

/*---------------------------------------------------------------------------*
 * ImrtTabuMemory
 *---------------------------------------------------------------------------*/

bool ImrtTabuMemory::tabu_check(emili::Solution *solution) {
  ImrtSolution *s = static_cast<ImrtSolution *>(solution);
  const std::vector<double> &x = s->getIntensities();
  for (int i = 0; i < count_; ++i) {
    int idx = (head_ - 1 - i + tenure_) % tenure_;
    if (memory_[idx] == x)
      return true;
  }
  return false;
}

void ImrtTabuMemory::forbid(emili::Solution *solution) {
  ImrtSolution *s = static_cast<ImrtSolution *>(solution);
  memory_[head_] = s->getIntensities();
  head_ = (head_ + 1) % tenure_;
  if (count_ < tenure_)
    ++count_;
}

} // namespace imrt
} // namespace emili
