#include "connections.h"

#include <cmath>
#include <limits>
#include <string>
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

  std::vector<std::pair<int, double>> get_connections(const int i) override;

 private:
  int n_dets = 0;

  int n_dets_prev = 0;

  int n_up = 0;

  int n_dn = 0;

  int cache_size = 0;

  std::vector<std::vector<std::pair<int, double>>> cached_connections;

  std::vector<uint8_t> cache_status;

  static constexpr uint8_t NOT_CACHED = 0;

  static constexpr uint8_t CACHED = 1;

  static constexpr uint8_t CACHE_OUTDATED = 2;

  static constexpr uint8_t CACHE_LIMIT_EXCEEDED = 4;

  // alpha and beta strings, O(n_dets).
  std::unordered_map<std::string, std::pair<std::vector<int>, std::vector<int>>>
      ab;

  // alpha-m1 and beta-m1 strings, O(n_dets * n_elecs).
  std::unordered_map<std::string, std::pair<std::vector<int>, std::vector<int>>>
      ab_m1;

  bool verbose = false;

  std::vector<std::vector<bool>> connected;

  std::vector<std::vector<bool>> one_up;

  Session* const session;

  AbstractSystem* const abstract_system;
};

ConnectionsImpl::ConnectionsImpl(
    Session* const session, AbstractSystem* const abstract_system)
    : session(session), abstract_system(abstract_system) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
  cache_size = config->get_int("cache_size");
  n_dets_prev = 0;
  connected.resize(omp_get_max_threads());
  one_up.resize(omp_get_max_threads());
}

void ConnectionsImpl::update() {
  // Alpha beta strings.
  const int n_dets_new = abstract_system->wf->terms_size();
  for (int i = n_dets; i < n_dets_new; i++) {
    const auto& det = abstract_system->wf->terms(i).det();
    ab[det.up().SerializeAsString()].first.push_back(i);
    ab[det.dn().SerializeAsString()].second.push_back(i);
  }
  if (verbose) printf("Number of unique alpha/beta: %'zu\n", ab.size());

  // Alpha beta minus one strings.
  for (int i = n_dets; i < n_dets_new; i++) {
    const auto& det = abstract_system->wf->terms(i).det();
    const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
    data::SpinDeterminant det_up(det.up());
    for (int j = 0; j < n_up; j++) {
      SpinDetUtil::set_occupation(&det_up, up_elecs[j], false);
      ab_m1[det_up.SerializeAsString()].first.push_back(i);
      SpinDetUtil::set_occupation(&det_up, up_elecs[j], true);
    }

    const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
    data::SpinDeterminant det_dn(det.dn());
    for (int j = 0; j < n_dn; j++) {
      SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], false);
      ab_m1[det_dn.SerializeAsString()].first.push_back(i);
      SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], true);
    }
  }
  if (verbose) printf("Number of unique alpha/beta-m1: %'zu\n", ab_m1.size());

#pragma omp parallel
  {
    // Connected and one up arrays.
    const int thread_id = omp_get_thread_num();
    connected[thread_id].assign(n_dets_new, false);
    one_up[thread_id].assign(n_dets_new, false);
  }

  // Cache.
  cached_connections.resize(n_dets_new);
  cache_status.resize(n_dets_new, NOT_CACHED);
  for (int i = 0; i < n_dets; i++) {
    if (cache_status[i] == CACHED) cache_status[i] = CACHE_OUTDATED;
  }

  n_dets_prev = n_dets;
  n_dets = n_dets_new;
};

std::vector<std::pair<int, double>> ConnectionsImpl::get_connections(
    const int i) {
  if (cache_status[i] == CACHED) {
    return cached_connections[i];
  }

  std::vector<std::pair<int, double>> res;
  if (cache_status[i] == CACHE_OUTDATED) {
    res = std::move(cached_connections[i]);
  }
  const int thread_id = omp_get_thread_num();
  const auto& det = abstract_system->wf->terms(i).det();

  // Two up/dn excitations.
  const auto& det_dn_code = det.dn().SerializeAsString();
  if (ab.find(det_dn_code) != ab.end()) {
    for (const int det_id : ab.find(det_dn_code)->second.second) {
      if (det_id < i) continue;
      if (cache_status[i] == CACHE_OUTDATED && det_id < n_dets_prev) continue;
      if (!connected[thread_id][det_id]) {
        const auto& det_id_det = abstract_system->wf->terms(det_id).det();
        const double H = abstract_system->hamiltonian(&det, &det_id_det);
        if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
        connected[thread_id][det_id] = true;
        res.push_back(std::make_pair(det_id, H));
      }
    }
  }
  const auto& det_up_code = det.up().SerializeAsString();
  if (ab.find(det_up_code) != ab.end()) {
    for (const int det_id : ab.find(det_up_code)->second.first) {
      if (det_id < i) continue;
      if (cache_status[i] == CACHE_OUTDATED && det_id < n_dets_prev) continue;
      if (!connected[thread_id][det_id]) {
        const auto& det_id_det = abstract_system->wf->terms(det_id).det();
        const double H = abstract_system->hamiltonian(&det, &det_id_det);
        if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
        connected[thread_id][det_id] = true;
        res.push_back(std::make_pair(det_id, H));
      }
    }
  }

  // One up one dn excitation.
  data::SpinDeterminant det_up(det.up());
  data::SpinDeterminant det_dn(det.dn());
  const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
  const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
  std::vector<int> one_ups;

  for (int k = 0; k < n_up; k++) {
    SpinDetUtil::set_occupation(&det_up, up_elecs[k], false);
    const auto& kv_up = ab_m1.find(det_up.SerializeAsString());
    if (kv_up != ab_m1.end()) {
      const auto& up_singles = kv_up->second.first;
      for (const int det_id : up_singles) {
        if (det_id < i) continue;
        if (cache_status[i] == CACHE_OUTDATED && det_id < n_dets_prev) continue;
        one_up[thread_id][det_id] = true;
        one_ups.push_back(det_id);
      }
    }
    SpinDetUtil::set_occupation(&det_up, up_elecs[k], true);
  }

  for (int k = 0; k < n_dn; k++) {
    SpinDetUtil::set_occupation(&det_dn, dn_elecs[k], false);
    const auto& kv_dn = ab_m1.find(det_dn.SerializeAsString());
    if (kv_dn != ab_m1.end()) {
      const auto& dn_singles = kv_dn->second.first;
      for (const int det_id : dn_singles) {
        if (det_id < i) continue;
        if (cache_status[i] == CACHE_OUTDATED && det_id < n_dets_prev) continue;
        if (one_up[thread_id][det_id] && !connected[thread_id][det_id]) {
          const auto& det_id_det = abstract_system->wf->terms(det_id).det();
          const double H = abstract_system->hamiltonian(&det, &det_id_det);
          if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
          connected[thread_id][det_id] = true;
          res.push_back(std::make_pair(det_id, H));
        }
        one_up[thread_id][det_id] = true;
        one_ups.push_back(det_id);
      }
    }
    SpinDetUtil::set_occupation(&det_dn, dn_elecs[k], true);
  }

  // Reset connected and return.
  for (const size_t det_id : one_ups) one_up[thread_id][det_id] = false;
  for (const auto& connection : res)
    connected[thread_id][connection.first] = false;

  const int n_connections = res.size();
  if (n_connections <= cache_size) {
    cached_connections[i] = res;
    cache_status[i] = CACHED;
  } else {
    cached_connections[i].clear();
    cache_status[i] = CACHE_LIMIT_EXCEEDED;
  }

  return res;
};

Connections* Injector::new_connections(
    Session* const session, AbstractSystem* const abstract_system) {
  return new ConnectionsImpl(session, abstract_system);
}
