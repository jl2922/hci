#include "heg_system.h"

#include <array>
#include <boost/functional/hash.hpp>
#include <climits>
#include <cmath>
#include <exception>
#include <unordered_map>
#include "../array_math.h"
#include "../injector.h"
#include "../solver/spin_det_util.h"
#include "k_points_util.h"
#include "omp.h"

class HEGSystemImpl : public HEGSystem {
 public:
  HEGSystemImpl(Session* const session);

  void setup(const double rcut) override;

  int get_n_orbitals(const double rcut) const override;

  double hamiltonian(
      const data::Determinant* const det_pq,
      const data::Determinant* const det_rs) const override;

  void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>&
          connected_det_handler) override;

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
      int,
      boost::hash<std::array<int8_t, 3>>>
      k_lut;

  std::unordered_map<
      std::array<int8_t, 3>,
      std::vector<std::pair<std::array<int8_t, 3>, double>>,
      boost::hash<std::array<int8_t, 3>>>
      same_spin_hci_queue;  // O(k_points^2).

  std::vector<std::pair<std::array<int8_t, 3>, double>>
      opposite_spin_hci_queue;  // O(k_points).

  std::vector<data::Determinant> tmp_dets;

  void setup_hci_queue();

  void evaluate_energy_hf();

  double hamiltonian_diagonal(const data::Determinant* const det) const;
};

HEGSystemImpl::HEGSystemImpl(Session* const session) : HEGSystem(session) {
  Config* const config = session->get_config();
  Parallel* const parallel = session->get_parallel();
  verbose = parallel->is_master();
  tmp_dets.resize(parallel->get_n_threads());
  r_s = config->get_double("r_s");
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
}

void HEGSystemImpl::setup(const double rcut) {
  assert(rcut > 0 && rcut < 32);
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
  for (int i = 0; i < n_k_points; i++) k_lut[k_points[i]] = i;

  // Setup HCI queue.
  setup_hci_queue();

  // Evaluate HF energy.
  evaluate_energy_hf();
}

int HEGSystemImpl::get_n_orbitals(const double rcut) const {
  return KPointsUtil::get_n_k_points(rcut) * 2;
}

double HEGSystemImpl::hamiltonian(
    const data::Determinant* const det_pq,
    const data::Determinant* const det_rs) const {
  if (det_pq == det_rs) {
    return hamiltonian_diagonal(det_pq);
  } else {
    const auto& eor_up = SpinDetUtil::get_eor(det_pq->up(), det_rs->up());
    const int n_eor_up = eor_up.size();
    if (n_eor_up > 4) return 0.0;
    const auto& eor_dn = SpinDetUtil::get_eor(det_pq->dn(), det_rs->dn());
    const int n_eor_dn = eor_dn.size();
    if (n_eor_up == 0 && n_eor_dn == 0) return hamiltonian_diagonal(det_pq);
    if (n_eor_up + n_eor_dn != 4) return 0.0;

    // Obtain p, q, s.
    bool k_p_set = false, k_r_set = false;
    int orb_p = 0, orb_r = 0, orb_s = 0;
    std::array<int8_t, 3> k_change;
    k_change.fill(0);
    for (const auto orb_i : eor_up) {
      if (SpinDetUtil::is_occupied(det_pq->up(), orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }
    for (const auto orb_i : eor_dn) {
      if (SpinDetUtil::is_occupied(det_pq->dn(), orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }

    // Check for momentum conservation.
    if (k_change != 0) return 0.0;

    double H = H_unit / squared_norm(k_points[orb_p] - k_points[orb_r]);
    if (n_eor_up != 2) {
      H -= H_unit / squared_norm(k_points[orb_p] - k_points[orb_s]);
    }

    const auto& get_gamma_exp = [&](const data::SpinDeterminant& spin_det,
                                    const std::vector<int>& eor) {
      int res = 0;
      for (const int orb : eor) {
        res += SpinDetUtil::get_n_lower_elecs(spin_det, orb);
      }
      return res;
    };

    const int gamma_exp = get_gamma_exp(det_pq->up(), eor_up) +
                          get_gamma_exp(det_pq->dn(), eor_dn) +
                          get_gamma_exp(det_rs->up(), eor_up) +
                          get_gamma_exp(det_rs->dn(), eor_dn);
    if ((gamma_exp & 1) == 1) H = -H;

    return H;
  }
};

double HEGSystemImpl::hamiltonian_diagonal(
    const data::Determinant* const det) const {
  double H = energy_hf;

  const auto& spin_correct = [&](const data::SpinDeterminant& spin_det) {
    const int n_hf_elecs = spin_det.n_hf_elecs();
    const int n_v_holes = spin_det.v_holes_size();
    const int n_c_elecs = spin_det.c_elecs_size();

    // One electron operators correction.
    for (const int p : spin_det.v_holes())
      H -= squared_norm(k_points[p] * k_unit) * 0.5;
    for (const int p : spin_det.c_elecs())
      H += squared_norm(k_points[p] * k_unit) * 0.5;

    // Two electron operators correction.
    std::vector<bool> hf_orbs(n_hf_elecs, true);
    for (const int p : spin_det.v_holes()) hf_orbs[p] = false;
    for (int p = 0; p < n_hf_elecs; p++) {
      if (!hf_orbs[p]) continue;
      const auto& k_p = k_points[p];
      for (const int q : spin_det.v_holes()) {
        const auto& k_q = k_points[q];
        H += H_unit / squared_norm(k_p - k_q);
      }
      for (const int q : spin_det.c_elecs()) {
        const auto& k_q = k_points[q];
        H -= H_unit / squared_norm(k_p - k_q);
      }
    }
    for (int ip = 0; ip < n_v_holes; ip++) {
      const auto& k_p = k_points[spin_det.v_holes(ip)];
      for (int iq = ip + 1; iq < n_v_holes; iq++) {
        const auto& k_q = k_points[spin_det.v_holes(iq)];
        H += H_unit / squared_norm(k_p - k_q);
      }
    }
    for (int ip = 0; ip < n_c_elecs; ip++) {
      const auto& k_p = k_points[spin_det.c_elecs(ip)];
      for (int iq = ip + 1; iq < n_c_elecs; iq++) {
        const auto& k_q = k_points[spin_det.c_elecs(iq)];
        H -= H_unit / squared_norm(k_p - k_q);
      }
    }
  };

  spin_correct(det->up());
  spin_correct(det->dn());

  return H;
}

void HEGSystemImpl::find_connected_dets(
    const data::Determinant* const det,
    const double eps,
    const std::function<void(const data::Determinant* const)>&
        connected_det_handler) {
  if (max_abs_H < eps) return;

  const int thread_id = omp_get_thread_num();
  auto& new_det = tmp_dets[thread_id];
  new_det = *det;

  const int dn_offset = k_points.size();

  const auto& is_occupied = [&](const data::Determinant* const target_det,
                                const int orb) {
    if (orb < dn_offset) {
      return SpinDetUtil::is_occupied(target_det->up(), orb);
    } else {
      return SpinDetUtil::is_occupied(target_det->dn(), orb - dn_offset);
    }
  };

  const auto& set_occupation =
      [&](data::Determinant* const target_det, const int orb, const bool occ) {
        if (orb < dn_offset) {
          SpinDetUtil::set_occupation(target_det->mutable_up(), orb, occ);
        } else {
          SpinDetUtil::set_occupation(
              target_det->mutable_dn(), orb - dn_offset, occ);
        }
      };

  const auto& pq_handler = [&](const int p, const int q) {
    int pp = p, qq = q;
    if (p >= dn_offset && q >= dn_offset) {
      pp -= dn_offset;
      qq -= dn_offset;
    } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
      pp = q - dn_offset;
      qq = p + dn_offset;
    }
    bool same_spin = false;
    const std::vector<std::pair<std::array<int8_t, 3>, double>>* items_ptr;
    if (pp < dn_offset && qq < dn_offset) {
      same_spin = true;
      const auto& diff_pq = k_points[qq] - k_points[pp];
      items_ptr = &(same_spin_hci_queue.find(diff_pq)->second);
    } else {
      items_ptr = &(opposite_spin_hci_queue);
    }
    const auto& items = *items_ptr;
    int qs_offset = 0;
    if (!same_spin) qs_offset = dn_offset;

    for (const auto& item : items) {
      if (item.second < eps) break;
      const auto& diff_pr = item.first;
      const auto& it_r = k_lut.find(diff_pr + k_points[pp]);
      if (it_r == k_lut.end()) continue;
      int r = it_r->second;
      const auto& it_s =
          k_lut.find(k_points[pp] + k_points[qq - qs_offset] - k_points[r]);
      if (it_s == k_lut.end()) continue;
      int s = it_s->second;
      if (same_spin && s < r) continue;
      s += qs_offset;
      if (p >= dn_offset && q >= dn_offset) {
        r += dn_offset;
        s += dn_offset;
      } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
        const int tmp = s;
        s = r + dn_offset;
        r = tmp - dn_offset;
      }

      if (!is_occupied(det, r) && !is_occupied(det, s)) {
        set_occupation(&new_det, p, false);
        set_occupation(&new_det, q, false);
        set_occupation(&new_det, r, true);
        set_occupation(&new_det, s, true);
        connected_det_handler(&new_det);
        set_occupation(&new_det, p, true);
        set_occupation(&new_det, q, true);
        set_occupation(&new_det, r, false);
        set_occupation(&new_det, s, false);
      }
    }
  };

  const auto& occ_up = SpinDetUtil::get_occupied_orbitals(det->up());
  const auto& occ_dn = SpinDetUtil::get_occupied_orbitals(det->dn());

  for (int i = 0; i < n_up; i++) {
    for (int j = i + 1; j < n_up; j++) {
      pq_handler(occ_up[i], occ_up[j]);
    }
  }
  for (int i = 0; i < n_dn; i++) {
    for (int j = i + 1; j < n_dn; j++) {
      pq_handler(occ_dn[i] + dn_offset, occ_dn[j] + dn_offset);
    }
  }
  for (int i = 0; i < n_up; i++) {
    for (int j = 0; j < n_dn; j++) {
      pq_handler(occ_up[i], occ_dn[j] + dn_offset);
    }
  }
};

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
