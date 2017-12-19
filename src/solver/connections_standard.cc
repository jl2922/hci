#include "connections.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "../injector.h"
#include "omp.h"
#include "spin_det_util.h"

class ConnectionsStandardImpl : public Connections {
 public:
  ConnectionsStandardImpl(
      Session* const session, AbstractSystem* const abstract_system);

  void update() override;

  void clear() override;

  std::vector<std::pair<int, double>> get_connections(
      const int det_id) override;

 private:
  int n_dets = 0;

  int n_dets_prev = 0;

  int n_up = 0;

  int n_dn = 0;

  std::vector<std::vector<std::pair<int, double>>> sparse_hamiltonian;

  bool verbose = false;

  std::unordered_map<std::string, int> alpha_to_id;

  std::unordered_map<std::string, int> beta_to_id;

  std::unordered_map<std::string, std::pair<std::vector<int>, std::vector<int>>>
      abm1_to_ab_ids;

  std::vector<std::vector<int>> alpha_id_to_single_ids;

  std::vector<std::vector<int>> beta_id_to_single_ids;

  // Sorted by unique beta id.
  std::vector<std::vector<int>> alpha_id_to_beta_ids;

  // Sorted by unique beta id.
  std::vector<std::vector<int>> alpha_id_to_det_ids;

  // Sorted by unique alpha id.
  std::vector<std::vector<int>> beta_id_to_alpha_ids;

  // Sorted by unique alpha id.
  std::vector<std::vector<int>> beta_id_to_det_ids;

  // Sort both vectors by the first vector.
  void sort_by_first(std::vector<int>& vec1, std::vector<int>& vec2);

  // Augment unique alphas/betas and alpha/beta to det info.
  void update_abdet();

  // Update unique alpha/beta minus one.
  void update_abm1();

  // Update alpha/beta singles lists.
  void update_absingles();

  void update_hamiltonian(const int det_id);
};

ConnectionsStandardImpl::ConnectionsStandardImpl(
    Session* const session, AbstractSystem* const abstract_system)
    : Connections(session, abstract_system) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
  n_dets = 0;
  n_dets_prev = 0;
}

void ConnectionsStandardImpl::clear() {
  n_dets = 0;
  n_dets_prev = 0;
  sparse_hamiltonian.clear();
  alpha_to_id.clear();
  beta_to_id.clear();
  abm1_to_ab_ids.clear();
  alpha_id_to_single_ids.clear();
  beta_id_to_single_ids.clear();
  alpha_id_to_beta_ids.clear();
  alpha_id_to_det_ids.clear();
  beta_id_to_alpha_ids.clear();
  beta_id_to_det_ids.clear();
}

void ConnectionsStandardImpl::update() {
  n_dets_prev = n_dets;
  n_dets = abstract_system->wf->terms_size();
  if (n_dets_prev == n_dets) return;

  // Construct helper lists.
  update_abdet();
  update_abm1();
  session->get_timer()->checkpoint("updated abm1");
  alpha_id_to_single_ids.clear();
  beta_id_to_single_ids.clear();
  update_absingles();
  session->get_timer()->checkpoint("updated absingles");
  abm1_to_ab_ids.clear();

  // Generate hamiltonian.
  sparse_hamiltonian.resize(n_dets);
  Parallel* const parallel = session->get_parallel();
  const int proc_id = parallel->get_proc_id();
  const int n_procs = parallel->get_n_procs();
  double target_progress = 0.01;
  if (verbose) {
    printf("Updating hamiltonian: ");
    fflush(stdout);
  }
#pragma omp parallel for schedule(dynamic, 10)
  for (int i = proc_id; i < n_dets; i += n_procs) {
    update_hamiltonian(i);
    if (verbose && omp_get_thread_num() == 0) {
      if (i > target_progress * n_dets) {
        printf("%.0f%% ", target_progress * 100);
        fflush(stdout);
        target_progress *= 2;
      }
    }
  }
  if (verbose) printf("\n");
  session->get_timer()->checkpoint("updated hamiltonian");
}

std::vector<std::pair<int, double>> ConnectionsStandardImpl::get_connections(
    const int det_id) {
  return sparse_hamiltonian[det_id];
}

void ConnectionsStandardImpl::update_abdet() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  for (int i = n_dets_prev; i < n_dets; i++) {
    const auto& det = abstract_system->wf->terms(i).det();

    // Obtain alpha id.
    const auto& alpha = det.up().SerializeAsString();
    int alpha_id;
    if (alpha_to_id.count(alpha) == 0) {
      alpha_id = alpha_to_id.size();
      alpha_to_id[alpha] = alpha_id;
      alpha_id_to_beta_ids.resize(alpha_id + 1);
      alpha_id_to_det_ids.resize(alpha_id + 1);
    } else {
      alpha_id = alpha_to_id[alpha];
    }

    // Obtain beta id.
    const auto& beta = det.dn().SerializeAsString();
    int beta_id;
    if (beta_to_id.count(beta) == 0) {
      beta_id = beta_to_id.size();
      beta_to_id[beta] = beta_id;
      beta_id_to_alpha_ids.resize(beta_id + 1);
      beta_id_to_det_ids.resize(beta_id + 1);
    } else {
      beta_id = beta_to_id[beta];
    }

    // Update alpha/beta to det info.
    alpha_id_to_beta_ids[alpha_id].push_back(beta_id);
    alpha_id_to_det_ids[alpha_id].push_back(i);
    beta_id_to_alpha_ids[beta_id].push_back(alpha_id);
    beta_id_to_det_ids[beta_id].push_back(i);
    updated_alphas.insert(alpha_id);
    updated_betas.insert(beta_id);
  }

  // Sort updated alpha/beta to det info.
  for (const int alpha_id : updated_alphas) {
    sort_by_first(
        alpha_id_to_beta_ids[alpha_id], alpha_id_to_det_ids[alpha_id]);
  }
  for (const int beta_id : updated_betas) {
    sort_by_first(beta_id_to_alpha_ids[beta_id], beta_id_to_det_ids[beta_id]);
  }
}

void ConnectionsStandardImpl::update_abm1() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  for (int i = n_dets_prev; i < n_dets; i++) {
    const auto& det = abstract_system->wf->terms(i).det();

    // Update alpha m1.
    const auto& alpha = det.up().SerializeAsString();
    const int alpha_id = alpha_to_id[alpha];
    if (updated_alphas.count(alpha_id) == 0) {
      const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
      data::SpinDeterminant det_up(det.up());
      for (int j = 0; j < n_up; j++) {
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], false);
        const auto& alpha_m1 = det_up.SerializeAsString();
        abm1_to_ab_ids[alpha_m1].first.push_back(alpha_id);
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], true);
      }
      updated_alphas.insert(alpha_id);
    }

    // Update beta m1.
    const auto& beta = det.dn().SerializeAsString();
    const int beta_id = beta_to_id[beta];
    if (updated_betas.count(beta_id) == 0) {
      const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
      data::SpinDeterminant det_dn(det.dn());
      for (int j = 0; j < n_dn; j++) {
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], false);
        const auto& beta_m1 = det_dn.SerializeAsString();
        abm1_to_ab_ids[beta_m1].second.push_back(beta_id);
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], true);
      }
      updated_betas.insert(beta_id);
    }
  }
  if (verbose) {
    printf("Outer size of abm1: %'zu\n", abm1_to_ab_ids.size());
    unsigned long long abm1_size = 0;
    for (const auto& item : abm1_to_ab_ids) {
      abm1_size += item.second.first.size();
      abm1_size += item.second.second.size();
    }
    printf("Full size of abm1: %'llu\n", abm1_size);
  }
}

void ConnectionsStandardImpl::update_absingles() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  alpha_id_to_single_ids.resize(alpha_to_id.size());
  beta_id_to_single_ids.resize(beta_to_id.size());
  for (int i = n_dets - 1; i >= 0; i--) {
    const auto& det = abstract_system->wf->terms(i).det();

    // Update alpha singles.
    const auto& alpha = det.up().SerializeAsString();
    const int alpha_id = alpha_to_id[alpha];
    if (updated_alphas.count(alpha_id) == 0) {
      const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
      data::SpinDeterminant det_up(det.up());
      for (int j = 0; j < n_up; j++) {
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], false);
        const auto& alpha_m1 = det_up.SerializeAsString();
        if (abm1_to_ab_ids.count(alpha_m1) == 1) {
          for (const int alpha_single : abm1_to_ab_ids[alpha_m1].first) {
            if (alpha_single == alpha_id) continue;
            if (i >= n_dets_prev && updated_alphas.count(alpha_single) == 1) {
              continue;
            }
            alpha_id_to_single_ids[alpha_id].push_back(alpha_single);
            alpha_id_to_single_ids[alpha_single].push_back(alpha_id);
          }
        }
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], true);
      }
      updated_alphas.insert(alpha_id);
    }

    // Update beta singles.
    const auto& beta = det.dn().SerializeAsString();
    const int beta_id = beta_to_id[beta];
    if (updated_betas.count(beta_id) == 0) {
      const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
      data::SpinDeterminant det_dn(det.dn());
      for (int j = 0; j < n_dn; j++) {
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], false);
        const auto& beta_m1 = det_dn.SerializeAsString();
        if (abm1_to_ab_ids.count(beta_m1) == 1) {
          for (const int beta_single : abm1_to_ab_ids[beta_m1].second) {
            if (beta_single == beta_id) continue;
            if (i >= n_dets_prev && updated_betas.count(beta_single) == 1) {
              continue;
            }
            beta_id_to_single_ids[beta_id].push_back(beta_single);
            beta_id_to_single_ids[beta_single].push_back(beta_id);
          }
        }
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], true);
      }
      updated_betas.insert(beta_id);
    }
  }

  // Sort updated alpha/beta singles and keep uniques.
  unsigned long long singles_cnt = 0;
  for (const int alpha_id : updated_alphas) {
    std::sort(
        alpha_id_to_single_ids[alpha_id].begin(),
        alpha_id_to_single_ids[alpha_id].end());
    singles_cnt += alpha_id_to_single_ids[alpha_id].size();
  }
  for (const int beta_id : updated_betas) {
    std::sort(
        beta_id_to_single_ids[beta_id].begin(),
        beta_id_to_single_ids[beta_id].end());
    singles_cnt += beta_id_to_single_ids[beta_id].size();
  }

  if (verbose) {
    printf(
        "Outer size of a/b singles: %'zu / %'zu\n",
        alpha_id_to_single_ids.size(),
        beta_id_to_single_ids.size());
    printf("Full size of absingles: %'llu\n", singles_cnt);
  }
}

void ConnectionsStandardImpl::update_hamiltonian(const int det_id) {
  std::vector<std::pair<int, double>>& res = sparse_hamiltonian[det_id];
  const auto& det = abstract_system->wf->terms(det_id).det();
  const bool is_new_det = det_id >= n_dets_prev;

  if (is_new_det) {
    const double H = abstract_system->hamiltonian(&det, &det);
    res.push_back(std::make_pair(det_id, H));
  }

  const int start_id = is_new_det ? det_id + 1 : n_dets_prev;

  // Single or double alpha excitations.
  const auto& beta = det.dn().SerializeAsString();
  const int beta_id = beta_to_id[beta];
  const auto& alpha_dets = beta_id_to_det_ids[beta_id];
  for (auto it = alpha_dets.begin(); it != alpha_dets.end(); it++) {
    const int alpha_det_id = *it;
    if (alpha_det_id < start_id) continue;
    const auto& alpha_det = abstract_system->wf->terms(alpha_det_id).det();
    const double H = abstract_system->hamiltonian(&det, &alpha_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    res.push_back(std::make_pair(alpha_det_id, H));
  }

  // Single or double beta excitations.
  const auto& alpha = det.up().SerializeAsString();
  const int alpha_id = alpha_to_id[alpha];
  const auto& beta_dets = alpha_id_to_det_ids[alpha_id];
  for (auto it = beta_dets.begin(); it != beta_dets.end(); it++) {
    const int beta_det_id = *it;
    if (beta_det_id < start_id) continue;
    const auto& beta_det = abstract_system->wf->terms(beta_det_id).det();
    const double H = abstract_system->hamiltonian(&det, &beta_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    res.push_back(std::make_pair(beta_det_id, H));
  }

  // Mixed double excitation.
  const auto& alpha_singles = alpha_id_to_single_ids[alpha_id];
  const auto& beta_singles = beta_id_to_single_ids[beta_id];
  for (const auto alpha_single : alpha_singles) {
    const auto& related_beta_ids = alpha_id_to_beta_ids[alpha_single];
    const auto& related_det_ids = alpha_id_to_det_ids[alpha_single];
    const int n_related_dets = related_beta_ids.size();
    int ptr = 0;
    for (auto it = beta_singles.begin(); it != beta_singles.end(); it++) {
      const int beta_single = *it;
      while (ptr < n_related_dets && related_beta_ids[ptr] < beta_single) {
        ptr++;
      }
      if (ptr == n_related_dets) break;

      if (related_beta_ids[ptr] == beta_single) {
        const int related_det_id = related_det_ids[ptr];
        ptr++;
        if (related_det_id < start_id) continue;
        const auto& related_det =
            abstract_system->wf->terms(related_det_id).det();
        const double H = abstract_system->hamiltonian(&det, &related_det);
        if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
        res.push_back(std::make_pair(related_det_id, H));
      }
    }
  }
}

void ConnectionsStandardImpl::sort_by_first(
    std::vector<int>& vec1, std::vector<int>& vec2) {
  std::vector<std::pair<int, int>> vec;
  const int n_vec = vec1.size();
  for (int i = 0; i < n_vec; i++) {
    vec.push_back(std::make_pair(vec1[i], vec2[i]));
  }
  std::sort(vec.begin(), vec.end(), [&](const auto& a, const auto& b) {
    return a.first < b.first;
  });
  for (int i = 0; i < n_vec; i++) {
    vec1[i] = vec[i].first;
    vec2[i] = vec[i].second;
  }
}

Connections* Injector::new_connections(
    Session* const session, AbstractSystem* const abstract_system) {
  return new ConnectionsStandardImpl(session, abstract_system);
}
