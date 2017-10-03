#include "heg_system.h"

#include <array>
#include <boost/functional/hash.hpp>
#include <climits>
#include <cmath>
#include <exception>
#include <unordered_map>
#include "../array_math.h"
#include "../injector.h"
#include "k_points_util.h"

class HEGSystemImpl : public HEGSystem {
 public:
  HEGSystemImpl(Session* const session);

  void setup(const double rcut) override;

  double hamiltonian(
      const data::Determinant* const det_pq,
      const data::Determinant* const det_rs) const override;

  void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>&
          connected_det_handler) const override;

 private:
  double rcut = 0.0;

  double r_s = 0.0;

  int n_up = 0;

  int n_dn = 0;

  double k_unit = 0.0;

  double H_unit = 0.0;

  double max_abs_H = 0.0;

  double energy_hf = 0.0;

  bool verbose = false;

  std::vector<std::array<int8_t, 3>> k_points;

  std::unordered_map<
      std::array<int8_t, 3>,
      size_t,
      boost::hash<std::array<int8_t, 3>>>
      k_lut;

  std::unordered_map<
      std::array<int8_t, 3>,
      std::vector<std::pair<std::array<int8_t, 3>, double>>,
      boost::hash<std::array<int8_t, 3>>>
      same_spin_hci_queue;  // O(k_points^2).

  std::vector<std::pair<std::array<int8_t, 3>, double>>
      opposite_spin_hci_queue;  // O(k_points).

  Session* const session;

  void setup_hci_queue();

  void evaluate_energy_hf();
};

HEGSystemImpl::HEGSystemImpl(Session* const session) : session(session) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  r_s = config->get_double("r_s");
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
}

void HEGSystemImpl::setup(const double rcut) {
  if (rcut <= 0 || rcut >= 127) {
    throw std::invalid_argument("rcut must be between 0 to 127");
  }
  this->rcut = rcut;

  constexpr double PI = 3.14159265358979323846;

  // Setup basic units.
  const double density = 3.0 / (4.0 * PI * std::pow(r_s, 3));
  const double cell_length = std::pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * PI / cell_length;
  H_unit = 1.0 / (PI * cell_length);

  // Setup K points and its look up table.
  k_points = KPointsUtil::generate_k_points(rcut);
  const int n_k_points = k_points.size();
  if (verbose) printf("Number of orbitals: %'d\n", n_k_points * 2);
  k_lut.clear();
  for (size_t i = 0; i < k_points.size(); i++) k_lut[k_points[i]] = i;

  // Setup HCI queue.
  setup_hci_queue();

  // Evaluate HF energy.
  evaluate_energy_hf();
}

double HEGSystemImpl::hamiltonian(
    const data::Determinant* const det_pq,
    const data::Determinant* const det_rs) const {
  return 0.0;
};

void HEGSystemImpl::find_connected_dets(
    const data::Determinant* const det,
    const double eps,
    const std::function<void(const data::Determinant* const)>&
        connected_det_handler) const {};

void HEGSystemImpl::setup_hci_queue() {
  same_spin_hci_queue.clear();
  opposite_spin_hci_queue.clear();
  max_abs_H = 0.0;

  // Common dependencies.
  const auto& k_diffs = KPointsUtil::get_k_diffs(k_points);
  const auto& sort_comparison = [](const auto& a, const auto& b) -> bool {
    return a.second > b.second;
  };

  // Same spin.
  const double diff_max_squared = 4.0 * rcut * rcut;
  for (const auto& diff_pq : k_diffs) {
    for (const auto& diff_pr : k_diffs) {
      const auto& diff_sr =
          diff_pr + diff_pr - diff_pq;  // Momentum conservation.
      if (diff_sr == 0 || squared_norm(diff_sr) > diff_max_squared) continue;
      const auto& diff_ps = diff_pr - diff_sr;
      if (diff_ps == 0) continue;
      if (std::abs(squared_norm(diff_pr) - squared_norm(diff_ps)) <
          std::numeric_limits<double>::epsilon()) {
        continue;
      }
      const double abs_H =
          std::abs(1.0 / squared_norm(diff_pr) - 1.0 / squared_norm(diff_ps));
      if (abs_H < std::numeric_limits<double>::epsilon()) continue;
      const auto& item = std::make_pair(diff_pr, abs_H * H_unit);
      same_spin_hci_queue[diff_pq].push_back(item);
    }
  }
  unsigned long long n_same_spin = 0;
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(items.begin(), items.end(), sort_comparison);
    max_abs_H = std::max(max_abs_H, items.front().second);
    n_same_spin += items.size();
  }
  if (verbose) {
    printf("Number of same spin hci queue items: %'llu\n", n_same_spin);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / squared_norm(diff_pr);
    if (abs_H < std::numeric_limits<double>::epsilon()) continue;
    const auto& item = std::make_pair(diff_pr, abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      sort_comparison);
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
  if (verbose) {
    unsigned long long n_opposite_spin = opposite_spin_hci_queue.size();
    printf("Number of opposite spin hci queue items: %'llu\n", n_opposite_spin);
  }
}

void HEGSystemImpl::evaluate_energy_hf() {
  double H = 0.0;

  // One electron operator.
  for (int p = 0; p < n_up; p++) H += squared_norm(k_points[p] * k_unit) * 0.5;
  for (int p = 0; p < n_dn; p++) H += squared_norm(k_points[p] * k_unit) * 0.5;

  // Two electrons operator.
  for (int p = 0; p < n_up; p++) {
    const auto& k_p = k_points[p];
    for (int q = p + 1; q < n_up; q++) {
      const auto& k_q = k_points[q];
      H -= H_unit / squared_norm(k_p - k_q);
    }
  }
  for (int p = 0; p < n_dn; p++) {
    const auto& k_p = k_points[p];
    for (int q = p + 1; q < n_dn; q++) {
      const auto& k_q = k_points[q];
      H -= H_unit / squared_norm(k_p - k_q);
    }
  }

  energy_hf = H;
  if (verbose) printf("HF energy: %.12f Ha\n", energy_hf);
}

HEGSystem* Injector::new_heg_system(Session* const session) {
  return new HEGSystemImpl(session);
}
