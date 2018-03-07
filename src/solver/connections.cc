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

class ConnectionsImpl : public Connections {
 public:
  ConnectionsImpl(
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

  std::vector<std::string> unique_alphas;

  std::vector<std::string> unique_betas;

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

  std::vector<data::Determinant> tmp_dets;

  std::vector<data::SpinDeterminant> tmp_spin_dets;

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

ConnectionsImpl::ConnectionsImpl(
    Session* const session, AbstractSystem* const abstract_system)
    : Connections(session, abstract_system) {
  Parallel* const parallel = session->get_parallel();
  Config* const config = session->get_config();
  verbose = parallel->is_master();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
  n_dets = 0;
  n_dets_prev = 0;
  const int n_threads = parallel->get_n_threads();
  tmp_dets.resize(n_threads * 2);
  tmp_spin_dets.resize(n_threads);
}

void ConnectionsImpl::clear() {
  n_dets = 0;
  n_dets_prev = 0;
  sparse_hamiltonian.clear();
  unique_alphas.clear();
  unique_betas.clear();
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

void ConnectionsImpl::update() {
  n_dets_prev = n_dets;
  n_dets = abstract_system->dets.size();
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

std::vector<std::pair<int, double>> ConnectionsImpl::get_connections(
    const int det_id) {
  return sparse_hamiltonian[det_id];
}

void ConnectionsImpl::update_abdet() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  auto& det = tmp_dets[0];
  for (int i = n_dets_prev; i < n_dets; i++) {
    det.ParseFromString(abstract_system->dets[i]);

    // Obtain alpha id.
    const auto& alpha = det.up().SerializeAsString();
    int alpha_id;
    if (alpha_to_id.count(alpha) == 0) {
      alpha_id = alpha_to_id.size();
      alpha_to_id[alpha] = alpha_id;
      unique_alphas.push_back(alpha);
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
      unique_betas.push_back(beta);
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

void ConnectionsImpl::update_abm1() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  auto& det = tmp_dets[0];
  auto& spin_det = tmp_spin_dets[0];
  for (int i = n_dets_prev; i < n_dets; i++) {
    det.ParseFromString(abstract_system->dets[i]);

    // Update alpha m1.
    const auto& alpha = det.up().SerializeAsString();
    const int alpha_id = alpha_to_id[alpha];
    if (updated_alphas.count(alpha_id) == 0) {
      const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
      spin_det = det.up();
      for (int j = 0; j < n_up; j++) {
        SpinDetUtil::set_occupation(&spin_det, up_elecs[j], false);
        const auto& alpha_m1 = spin_det.SerializeAsString();
        abm1_to_ab_ids[alpha_m1].first.push_back(alpha_id);
        SpinDetUtil::set_occupation(&spin_det, up_elecs[j], true);
      }
      updated_alphas.insert(alpha_id);
    }

    // Update beta m1.
    const auto& beta = det.dn().SerializeAsString();
    const int beta_id = beta_to_id[beta];
    if (updated_betas.count(beta_id) == 0) {
      const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
      spin_det = det.dn();
      for (int j = 0; j < n_dn; j++) {
        SpinDetUtil::set_occupation(&spin_det, dn_elecs[j], false);
        const auto& beta_m1 = spin_det.SerializeAsString();
        abm1_to_ab_ids[beta_m1].second.push_back(beta_id);
        SpinDetUtil::set_occupation(&spin_det, dn_elecs[j], true);
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

void ConnectionsImpl::update_absingles() {
  std::unordered_set<int> updated_alphas;
  std::unordered_set<int> updated_betas;
  alpha_id_to_single_ids.resize(alpha_to_id.size());
  beta_id_to_single_ids.resize(beta_to_id.size());
  auto& det = tmp_dets[0];

  for (int i = n_dets_prev; i < n_dets; i++) {
    det.ParseFromString(abstract_system->dets[i]);

    const auto& alpha = det.up().SerializeAsString();
    const int alpha_id = alpha_to_id[alpha];
    updated_alphas.insert(alpha_id);

    const auto& beta = det.dn().SerializeAsString();
    const int beta_id = beta_to_id[beta];
    updated_betas.insert(beta_id);
  }

  const int n_unique_alphas = alpha_to_id.size();
  const int n_unique_betas = beta_to_id.size();

  std::vector<omp_lock_t> locks;
  const int n_locks = std::max(n_unique_alphas, n_unique_betas);
  locks.resize(n_locks);
  for (auto& lock : locks) omp_init_lock(&lock);

#pragma omp parallel for schedule(static, 1)
  for (int alpha_id = 0; alpha_id < n_unique_alphas; alpha_id++) {
    const auto& alpha = unique_alphas[alpha_id];
    const int thread_id = omp_get_thread_num();
    auto& spin_det = tmp_spin_dets[thread_id];
    spin_det.ParseFromString(alpha);
    const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(spin_det);
    for (int j = 0; j < n_up; j++) {
      SpinDetUtil::set_occupation(&spin_det, up_elecs[j], false);
      const auto& alpha_m1 = spin_det.SerializeAsString();
      if (abm1_to_ab_ids.count(alpha_m1) == 1) {
        for (const int alpha_single : abm1_to_ab_ids[alpha_m1].first) {
          if (alpha_single == alpha_id) continue;
          if (alpha_id > alpha_single && updated_alphas.count(alpha_id) &&
              updated_alphas.count(alpha_single)) {
            continue;
          }
          omp_set_lock(&locks[alpha_id]);
          alpha_id_to_single_ids[alpha_id].push_back(alpha_single);
          omp_unset_lock(&locks[alpha_id]);
          omp_set_lock(&locks[alpha_single]);
          alpha_id_to_single_ids[alpha_single].push_back(alpha_id);
          omp_unset_lock(&locks[alpha_single]);
        }
      }
      SpinDetUtil::set_occupation(&spin_det, up_elecs[j], true);
    }
  }

#pragma omp parallel for schedule(static, 1)
  for (int beta_id = 0; beta_id < n_unique_betas; beta_id++) {
    const auto& beta = unique_betas[beta_id];
    const int thread_id = omp_get_thread_num();
    auto& spin_det = tmp_spin_dets[thread_id];
    spin_det.ParseFromString(beta);
    const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(spin_det);
    for (int j = 0; j < n_dn; j++) {
      SpinDetUtil::set_occupation(&spin_det, dn_elecs[j], false);
      const auto& beta_m1 = spin_det.SerializeAsString();
      if (abm1_to_ab_ids.count(beta_m1) == 1) {
        for (const int beta_single : abm1_to_ab_ids[beta_m1].second) {
          if (beta_single == beta_id) continue;
          if (beta_id > beta_single && updated_betas.count(beta_id) &&
              updated_betas.count(beta_single)) {
            continue;
          }
          omp_set_lock(&locks[beta_id]);
          beta_id_to_single_ids[beta_id].push_back(beta_single);
          omp_unset_lock(&locks[beta_id]);
          omp_set_lock(&locks[beta_single]);
          beta_id_to_single_ids[beta_single].push_back(beta_id);
          omp_unset_lock(&locks[beta_single]);
        }
      }
      SpinDetUtil::set_occupation(&spin_det, dn_elecs[j], true);
    }
  }

  for (auto& lock : locks) omp_destroy_lock(&lock);

  // Sort updated alpha/beta singles and keep uniques.
  unsigned long long singles_cnt = 0;
#pragma omp parallel for schedule(static, 1) reduction(+ : singles_cnt)
  for (int alpha_id = 0; alpha_id < n_unique_alphas; alpha_id++) {
    std::sort(
        alpha_id_to_single_ids[alpha_id].begin(),
        alpha_id_to_single_ids[alpha_id].end());
    singles_cnt += alpha_id_to_single_ids[alpha_id].size();
  }
#pragma omp parallel for schedule(static, 1) reduction(+ : singles_cnt)
  for (int beta_id = 0; beta_id < n_unique_betas; beta_id++) {
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

void ConnectionsImpl::update_hamiltonian(const int det_id) {
  std::vector<std::pair<int, double>>& res = sparse_hamiltonian[det_id];
  const int thread_id = omp_get_thread_num();
  auto& det = tmp_dets[thread_id * 2];
  const double coef = abstract_system->coefs[det_id];
  auto& connected_det = tmp_dets[thread_id * 2 + 1];
  det.ParseFromString(abstract_system->dets[det_id]);
  const bool is_new_det = det_id >= n_dets_prev;

  if (is_new_det) {
    const double H = abstract_system->hamiltonian(&det, &det);
    res.push_back(std::make_pair(det_id, H));
  }
  // if (det_id >= 57315) return;

  const int start_id = is_new_det ? det_id + 1 : n_dets_prev;

  // Single or double alpha excitations.
  const auto& beta = det.dn().SerializeAsString();
  const int beta_id = beta_to_id[beta];
  const auto& alpha_dets = beta_id_to_det_ids[beta_id];
  for (auto it = alpha_dets.begin(); it != alpha_dets.end(); it++) {
    const int alpha_det_id = *it;
    // if (alpha_det_id >= 57315) continue;
    if (alpha_det_id < start_id) continue;
    connected_det.ParseFromString(abstract_system->dets[alpha_det_id]);
    const double H = abstract_system->hamiltonian(&det, &connected_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    if (det_id >= 57315 && alpha_det_id < 57315 && std::abs(abstract_system->coefs[alpha_det_id] * H) < 0.00002) continue;
    if (det_id < 57315 && alpha_det_id >= 57315 && std::abs(coef * H) < 0.00002) continue;
    if (det_id >= 57315 && alpha_det_id >= 57315) continue;
    if (det_id >= 57315 && alpha_det_id >= 57315 && std::abs(coef * H) < 0.00002 && 
         std::abs(abstract_system->coefs[alpha_det_id] * H) < 0.00002) continue;
    res.push_back(std::make_pair(alpha_det_id, H));
  }

  // Single or double beta excitations.
  const auto& alpha = det.up().SerializeAsString();
  const int alpha_id = alpha_to_id[alpha];
  const auto& beta_dets = alpha_id_to_det_ids[alpha_id];
  for (auto it = beta_dets.begin(); it != beta_dets.end(); it++) {
    const int beta_det_id = *it;
    // if (beta_det_id >= 57315) continue;
    if (beta_det_id < start_id) continue;
    connected_det.ParseFromString(abstract_system->dets[beta_det_id]);
    const double H = abstract_system->hamiltonian(&det, &connected_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    if (det_id >= 57315 && beta_det_id < 57315 && std::abs(abstract_system->coefs[beta_det_id] * H) < 0.00002) continue;
    if (det_id < 57315 && beta_det_id >= 57315 && std::abs(coef * H) < 0.00002) continue;
    if (det_id >= 57315 && beta_det_id >= 57315) continue;
    if (det_id >= 57315 && beta_det_id >= 57315 && std::abs(coef * H) < 0.00002 && 
         std::abs(abstract_system->coefs[beta_det_id] * H) < 0.00002) continue;
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
        // if (related_det_id >= 57315) continue;
        if (related_det_id < start_id) continue;
        connected_det.ParseFromString(abstract_system->dets[related_det_id]);
        const double H = abstract_system->hamiltonian(&det, &connected_det);
        if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
        if (det_id >= 57315 && related_det_id < 57315 && std::abs(abstract_system->coefs[related_det_id] * H) < 0.00002) continue;
        if (det_id < 57315 && related_det_id >= 57315 && std::abs(coef * H) < 0.00002) continue;
        if (det_id >= 57315 && related_det_id >= 57315) continue;
        if (det_id >= 57315 && related_det_id >= 57315 && std::abs(coef * H) < 0.00002 && 
             std::abs(abstract_system->coefs[related_det_id] * H) < 0.00002) continue;
        res.push_back(std::make_pair(related_det_id, H));
      }
    }
  }
}

void ConnectionsImpl::sort_by_first(
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
  return new ConnectionsImpl(session, abstract_system);
}
