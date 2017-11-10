#include "connections.h"

#include <algorithm>
#include <boost/functional/hash/hash.hpp>
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

class ConnectionsABSinglesImpl : public Connections {
 public:
  ConnectionsABSinglesImpl(
      Session* const session, AbstractSystem* const abstract_system);

  void update() override;

  void clear() override;

  std::vector<std::pair<int, double>> get_connections(const int i) override;

 private:
  int n_dets = 0;

  int n_dets_prev = 0;

  int n_up = 0;

  int n_dn = 0;

  int cache_size = 0;

  std::vector<std::vector<std::pair<int, double>>> cached_connections;

  static constexpr uint8_t NOT_CACHED = 0;

  static constexpr uint8_t CACHED = 1;

  static constexpr uint8_t CACHE_OUTDATED = 2;

  static constexpr uint8_t CACHE_LIMIT_EXCEEDED = 3;

  std::vector<uint8_t> cache_status;

  bool verbose = false;

  std::unordered_map<std::string, int> unique_alphas;

  std::unordered_map<std::string, int> unique_betas;

  std::unordered_map<std::string, std::pair<std::vector<int>, std::vector<int>>>
      unique_ab_m1;

  std::unordered_map<int, std::vector<int>> singles_from_alpha;

  std::unordered_map<int, std::vector<int>> singles_from_beta;

  std::unordered_map<int, std::vector<int>> alpha_major_to_det;

  std::unordered_map<int, std::vector<int>> beta_major_to_det;

  std::vector<std::vector<bool>> one_up;
};

constexpr uint8_t ConnectionsABSinglesImpl::NOT_CACHED;
constexpr uint8_t ConnectionsABSinglesImpl::CACHED;
constexpr uint8_t ConnectionsABSinglesImpl::CACHE_OUTDATED;
constexpr uint8_t ConnectionsABSinglesImpl::CACHE_LIMIT_EXCEEDED;

ConnectionsABSinglesImpl::ConnectionsABSinglesImpl(
    Session* const session, AbstractSystem* const abstract_system)
    : Connections(session, abstract_system) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
  cache_size = config->get_int("cache_size");
  n_dets = 0;
  n_dets_prev = 0;
  one_up.resize(omp_get_max_threads());
}

void ConnectionsABSinglesImpl::clear() {
  n_dets = 0;
  n_dets_prev = 0;
  cache_status.clear();
  cached_connections.clear();
  unique_alphas.clear();
  unique_betas.clear();
  unique_ab_m1.clear();
  singles_from_alpha.clear();
  singles_from_beta.clear();
  alpha_major_to_det.clear();
  beta_major_to_det.clear();
  one_up.clear();
}

void ConnectionsABSinglesImpl::update() {
  n_dets_prev = n_dets;
  n_dets = abstract_system->wf->terms_size();
  if (n_dets_prev == n_dets) return;

  for (int i = n_dets_prev; i < n_dets; i++) {
    const auto& det = abstract_system->wf->terms(i).det();

    const auto& alpha = det.up().SerializeAsString();
    const auto& beta = det.dn().SerializeAsString();
    int alpha_id = -1;
    int beta_id = -1;

    // Process alpha.
    if (unique_alphas.count(alpha) == 1) {
      alpha_id = unique_alphas[alpha];
    } else {
      alpha_id = unique_alphas.size();
      unique_alphas[alpha] = alpha_id;

      // Setup connections.
      const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
      data::SpinDeterminant det_up(det.up());
      std::unordered_set<int> alpha_singles_set;
      for (int j = 0; j < n_up; j++) {
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], false);
        const auto& alpha_m1 = det_up.SerializeAsString();
        for (const int alpha_single : unique_ab_m1[alpha_m1].first) {
          if (alpha_singles_set.count(alpha_single) == 0) {
            singles_from_alpha[alpha_id].push_back(alpha_single);
            alpha_singles_set.insert(alpha_single);
          }
          if (singles_from_alpha[alpha_single].empty() ||
              singles_from_alpha[alpha_single].back() != alpha_id) {
            singles_from_alpha[alpha_single].push_back(alpha_id);
          }
        }
        unique_ab_m1[alpha_m1].first.push_back(alpha_id);
        SpinDetUtil::set_occupation(&det_up, up_elecs[j], true);
      }
    }

    // Process beta.
    if (unique_betas.count(beta) == 1) {
      beta_id = unique_betas[beta];
    } else {
      beta_id = unique_betas.size();
      unique_betas[beta] = beta_id;

      // Setup connections.
      const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
      data::SpinDeterminant det_dn(det.dn());
      std::unordered_set<int> single_betas_set;
      for (int j = 0; j < n_dn; j++) {
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], false);
        const auto& beta_m1 = det_dn.SerializeAsString();
        for (const int beta_single : unique_ab_m1[beta_m1].second) {
          if (single_betas_set.count(beta_single) == 0) {
            singles_from_beta[beta_id].push_back(beta_single);
            single_betas_set.insert(beta_single);
          }
          if (singles_from_beta[beta_single].empty() ||
              singles_from_beta[beta_single].back() != beta_id) {
            singles_from_beta[beta_single].push_back(beta_id);
          }
        }
        unique_ab_m1[beta_m1].second.push_back(beta_id);
        SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], true);
      }
    }

    alpha_major_to_det[alpha_id].push_back(i);
    beta_major_to_det[beta_id].push_back(i);
  }

#pragma omp parallel
  {
    // Connected and one up arrays.
    const int thread_id = omp_get_thread_num();
    one_up[thread_id].assign(n_dets, false);
  }

  // Cache.
  cached_connections.resize(n_dets);
  cache_status.resize(n_dets, NOT_CACHED);
  for (int i = 0; i < n_dets_prev; i++) {
    if (cache_status[i] == CACHED) cache_status[i] = CACHE_OUTDATED;
  }
}

std::vector<std::pair<int, double>> ConnectionsABSinglesImpl::get_connections(
    const int i) {
  if (cache_status[i] == CACHED) {
    return cached_connections[i];
  }

  std::vector<std::pair<int, double>> res;

  // If already cached in the previous iteration, start from previous results.
  if (cache_status[i] == CACHE_OUTDATED) {
    res = std::move(cached_connections[i]);
  }

  const auto& det = abstract_system->wf->terms(i).det();
  const int thread_id = omp_get_thread_num();

  // Diagonal.
  if (cache_status[i] != CACHE_OUTDATED) {
    // Otherwise, diagonal term must be in the cache already.
    const double H = abstract_system->hamiltonian(&det, &det);
    res.push_back(std::make_pair(i, H));
  }

  const int start_id = cache_status[i] == CACHE_OUTDATED ? n_dets_prev : i + 1;

  // Two up excitation.
  const auto& beta = det.dn().SerializeAsString();
  const int beta_id = unique_betas[beta];
  const auto& two_ups = beta_major_to_det[beta_id];
  const auto& start_it_two_ups =
      std::lower_bound(two_ups.begin(), two_ups.end(), start_id);
  for (auto it = start_it_two_ups; it != two_ups.end(); it++) {
    const int det_id = *it;
    const auto& det_id_det = abstract_system->wf->terms(det_id).det();
    const double H = abstract_system->hamiltonian(&det, &det_id_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    res.push_back(std::make_pair(det_id, H));
  }

  // Two down exciation.
  const auto& alpha = det.up().SerializeAsString();
  const int alpha_id = unique_alphas[alpha];
  const auto& two_dns = alpha_major_to_det[alpha_id];
  const auto& start_it_two_dns =
      std::lower_bound(two_dns.begin(), two_dns.end(), start_id);
  for (auto it = start_it_two_dns; it != two_dns.end(); it++) {
    const int det_id = *it;
    const auto& det_id_det = abstract_system->wf->terms(det_id).det();
    const double H = abstract_system->hamiltonian(&det, &det_id_det);
    if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
    res.push_back(std::make_pair(det_id, H));
  }

  // One up one down exciation.
  std::vector<int> one_ups;
  const auto& single_alphas = singles_from_alpha[alpha_id];
  for (auto it = single_alphas.begin(); it != single_alphas.end(); it++) {
    const int single_alpha_id = *it;
    const auto& alpha_dets = alpha_major_to_det[single_alpha_id];

    const auto& start_it =
        std::lower_bound(alpha_dets.begin(), alpha_dets.end(), start_id);
    for (auto it_det = start_it; it_det < alpha_dets.end(); it_det++) {
      const int det_id = *it_det;
      one_up[thread_id][det_id] = true;
      one_ups.push_back(det_id);
    }
  }
  const auto& single_betas = singles_from_beta[beta_id];
  for (auto it = single_betas.begin(); it != single_betas.end(); it++) {
    const int single_beta_id = *it;
    const auto& beta_dets = beta_major_to_det[single_beta_id];

    const auto& start_it =
        std::lower_bound(beta_dets.begin(), beta_dets.end(), start_id);
    for (auto it_det = start_it; it_det < beta_dets.end(); it_det++) {
      const int det_id = *it_det;
      if (one_up[thread_id][det_id]) {
        const auto& det_id_det = abstract_system->wf->terms(det_id).det();
        const double H = abstract_system->hamiltonian(&det, &det_id_det);
        if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
        res.push_back(std::make_pair(det_id, H));
      }
    }
  }
  for (const size_t det_id : one_ups) one_up[thread_id][det_id] = false;

  // Cache if within threshold.
  const int n_connections = res.size();
  if (n_connections <= cache_size) {
    cached_connections[i] = res;
    cache_status[i] = CACHED;
  } else {
    cached_connections[i].clear();
    cache_status[i] = CACHE_LIMIT_EXCEEDED;
  }

  return res;
}

// Connections* Injector::new_connections(
//     Session* const session, AbstractSystem* const abstract_system) {
//   return new ConnectionsABSinglesImpl(session, abstract_system);
// }