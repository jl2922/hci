// #include "connections.h"

// #include <algorithm>
// #include <boost/functional/hash/hash.hpp>
// #include <cmath>
// #include <limits>
// #include <string>
// #include <unordered_map>
// #include <unordered_set>
// #include <utility>
// #include <vector>
// #include "../injector.h"
// #include "omp.h"
// #include "spin_det_util.h"

// class ConnectionsSandeepImpl : public Connections {
//  public:
//   ConnectionsSandeepImpl(
//       Session* const session, AbstractSystem* const abstract_system);

//   void update() override;

//   void clear() override;

//   std::vector<std::pair<int, double>> get_connections(const int i) override;

//  private:
//   int n_dets = 0;

//   int n_dets_prev = 0;

//   int n_up = 0;

//   int n_dn = 0;

//   // Maximum number of connetions to cache for each det.
//   size_t cache_size = 0;

//   // Maximum number of single connections to cache for each unique a/b.
//   size_t singles_cache_size = 0;

//   std::vector<std::vector<std::pair<int, double>>> cached_connections;

//   static constexpr uint8_t NOT_CACHED = 0;

//   static constexpr uint8_t CACHED = 1;

//   static constexpr uint8_t CACHE_OUTDATED = 2;

//   static constexpr uint8_t CACHE_OVERFLOW = 3;

//   std::vector<uint8_t> cache_status;

//   bool verbose = false;

//   std::unordered_map<std::string, int> unique_alphas;

//   std::unordered_map<std::string, int> unique_betas;

//   std::unordered_map<std::string, std::pair<std::vector<int>,
//   std::vector<int>>>
//       unique_ab_m1;

//   std::unordered_map<int, std::vector<int>> singles_from_alpha;

//   std::unordered_map<int, std::vector<int>> singles_from_beta;

//   std::vector<uint8_t> singles_alpha_cache_status;

//   std::vector<uint8_t> singles_beta_cache_status;

//   // Sorted by unique beta id.
//   std::unordered_map<int, std::vector<int>> alpha_major_to_beta;

//   // Sorted by unique beta id.
//   std::unordered_map<int, std::vector<int>> alpha_major_to_det;

//   // Sorted by det id.
//   std::unordered_map<int, std::vector<int>> alpha_major_to_det2;

//   // Sorted by unique alpha id.
//   std::unordered_map<int, std::vector<int>> beta_major_to_alpha;

//   // Sorted by unique alpha id.
//   std::unordered_map<int, std::vector<int>> beta_major_to_det;

//   // Sorted by det id.
//   std::unordered_map<int, std::vector<int>> beta_major_to_det2;

//   // Reusing an array of n_dets false values for efficiency.
//   std::vector<std::vector<bool>> one_up;

//   // Sort both vectors by the first vector.
//   void sort_by_first(std::vector<int>& vec1, std::vector<int>& vec2);

//   // Get singles a/b from cache or calculate on the fly if not cached..
//   std::vector<int> get_alpha_singles(
//       const int i, const data::SpinDeterminant& det_up) const;

//   std::vector<int> get_beta_singles(
//       const int i, const data::SpinDeterminant& det_dn) const;
// };

// constexpr uint8_t ConnectionsSandeepImpl::NOT_CACHED;
// constexpr uint8_t ConnectionsSandeepImpl::CACHED;
// constexpr uint8_t ConnectionsSandeepImpl::CACHE_OUTDATED;
// constexpr uint8_t ConnectionsSandeepImpl::CACHE_OVERFLOW;

// ConnectionsSandeepImpl::ConnectionsSandeepImpl(
//     Session* const session, AbstractSystem* const abstract_system)
//     : Connections(session, abstract_system) {
//   verbose = session->get_parallel()->is_master();
//   Config* const config = session->get_config();
//   n_up = config->get_int("n_up");
//   n_dn = config->get_int("n_dn");
//   cache_size = config->get_int("cache_size");
//   singles_cache_size = config->get_int("singles_cache_size");
//   n_dets = 0;
//   n_dets_prev = 0;
//   one_up.resize(omp_get_max_threads());
// }

// void ConnectionsSandeepImpl::clear() {
//   n_dets = 0;
//   n_dets_prev = 0;
//   cache_status.clear();
//   cached_connections.clear();
//   unique_alphas.clear();
//   unique_betas.clear();
//   unique_ab_m1.clear();
//   singles_from_alpha.clear();
//   singles_from_beta.clear();
//   singles_alpha_cache_status.clear();
//   singles_beta_cache_status.clear();
//   alpha_major_to_beta.clear();
//   alpha_major_to_det.clear();
//   alpha_major_to_det2.clear();
//   beta_major_to_alpha.clear();
//   beta_major_to_det.clear();
//   beta_major_to_det2.clear();
// }

// void ConnectionsSandeepImpl::update() {
//   n_dets_prev = n_dets;
//   n_dets = abstract_system->wf->terms_size();
//   if (n_dets_prev == n_dets) return;

//   singles_alpha_cache_status.resize(n_dets, CACHED);
//   singles_beta_cache_status.resize(n_dets, CACHED);

//   std::unordered_set<int> changed_alphas;
//   std::unordered_set<int> changed_betas;

//   for (int i = n_dets_prev; i < n_dets; i++) {
//     const auto& det = abstract_system->wf->terms(i).det();

//     const auto& alpha = det.up().SerializeAsString();
//     const auto& beta = det.dn().SerializeAsString();
//     int alpha_id = -1;
//     int beta_id = -1;

//     // Process alpha.
//     if (unique_alphas.count(alpha) == 1) {
//       alpha_id = unique_alphas[alpha];
//     } else {
//       alpha_id = unique_alphas.size();
//       unique_alphas[alpha] = alpha_id;

//       // Setup connections.
//       const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det.up());
//       data::SpinDeterminant det_up(det.up());
//       std::unordered_set<int> alpha_singles_set;
//       for (int j = 0; j < n_up; j++) {
//         SpinDetUtil::set_occupation(&det_up, up_elecs[j], false);
//         const auto& alpha_m1 = det_up.SerializeAsString();
//         for (const int alpha_single : unique_ab_m1[alpha_m1].first) {
//           if (singles_alpha_cache_status[alpha_id] != CACHE_OVERFLOW &&
//               alpha_singles_set.count(alpha_single) == 0) {
//             singles_from_alpha[alpha_id].push_back(alpha_single);
//             alpha_singles_set.insert(alpha_single);
//             if (singles_from_alpha[alpha_id].size() > singles_cache_size) {
//               singles_alpha_cache_status[alpha_id] = CACHE_OVERFLOW;
//               singles_from_alpha[alpha_id].clear();
//             }
//           }
//           if (singles_alpha_cache_status[alpha_single] != CACHE_OVERFLOW &&
//               (singles_from_alpha[alpha_single].empty() ||
//                singles_from_alpha[alpha_single].back() != alpha_id)) {
//             singles_from_alpha[alpha_single].push_back(alpha_id);
//             if (singles_from_alpha[alpha_single].size() > singles_cache_size)
//             {
//               singles_alpha_cache_status[alpha_single] = CACHE_OVERFLOW;
//               singles_from_alpha[alpha_single].clear();
//             }
//           }
//         }
//         unique_ab_m1[alpha_m1].first.push_back(alpha_id);
//         SpinDetUtil::set_occupation(&det_up, up_elecs[j], true);
//       }
//       std::sort(
//           singles_from_alpha[alpha_id].begin(),
//           singles_from_alpha[alpha_id].end());
//     }

//     // Process beta.
//     if (unique_betas.count(beta) == 1) {
//       beta_id = unique_betas[beta];
//     } else {
//       beta_id = unique_betas.size();
//       unique_betas[beta] = beta_id;

//       // Setup connections.
//       const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det.dn());
//       data::SpinDeterminant det_dn(det.dn());
//       std::unordered_set<int> single_betas_set;
//       for (int j = 0; j < n_dn; j++) {
//         SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], false);
//         const auto& beta_m1 = det_dn.SerializeAsString();
//         for (const int beta_single : unique_ab_m1[beta_m1].second) {
//           if (singles_beta_cache_status[beta_id] != CACHE_OVERFLOW &&
//               single_betas_set.count(beta_single) == 0) {
//             singles_from_beta[beta_id].push_back(beta_single);
//             single_betas_set.insert(beta_single);
//             if (singles_from_beta[beta_id].size() > singles_cache_size) {
//               singles_beta_cache_status[beta_id] = CACHE_OVERFLOW;
//               singles_from_beta[beta_id].clear();
//             }
//           }
//           if (singles_beta_cache_status[beta_single] != CACHE_OVERFLOW &&
//               (singles_from_beta[beta_single].empty() ||
//                singles_from_beta[beta_single].back() != beta_id)) {
//             singles_from_beta[beta_single].push_back(beta_id);
//             if (singles_from_beta[beta_single].size() > singles_cache_size) {
//               singles_beta_cache_status[beta_single] = CACHE_OVERFLOW;
//               singles_from_beta[beta_single].clear();
//             }
//           }
//         }
//         unique_ab_m1[beta_m1].second.push_back(beta_id);
//         SpinDetUtil::set_occupation(&det_dn, dn_elecs[j], true);
//       }
//       std::sort(
//           singles_from_beta[beta_id].begin(),
//           singles_from_beta[beta_id].end());
//     }

//     alpha_major_to_beta[alpha_id].push_back(beta_id);
//     alpha_major_to_det[alpha_id].push_back(i);
//     alpha_major_to_det2[alpha_id].push_back(i);
//     beta_major_to_alpha[beta_id].push_back(alpha_id);
//     beta_major_to_det[beta_id].push_back(i);
//     beta_major_to_det2[beta_id].push_back(i);

//     if (changed_alphas.count(alpha_id) == 0) {
//       changed_alphas.insert(alpha_id);
//     }

//     if (changed_betas.count(beta_id) == 0) {
//       changed_betas.insert(beta_id);
//     }
//   }

//   for (const int alpha_id : changed_alphas) {
//     sort_by_first(alpha_major_to_beta[alpha_id],
//     alpha_major_to_det[alpha_id]);
//   }
//   for (const int beta_id : changed_betas) {
//     sort_by_first(beta_major_to_alpha[beta_id], beta_major_to_det[beta_id]);
//   }

//   // Cache.
//   cached_connections.resize(n_dets);
//   cache_status.resize(n_dets, NOT_CACHED);
//   for (int i = 0; i < n_dets_prev; i++) {
//     if (cache_status[i] == CACHED) cache_status[i] = CACHE_OUTDATED;
//   }

// #pragma omp parallel
//   {
//     // Connected and one up arrays.
//     const int thread_id = omp_get_thread_num();
//     one_up[thread_id].assign(n_dets, false);
//   }
// }

// std::vector<int> ConnectionsSandeepImpl::get_alpha_singles(
//     const int i, const data::SpinDeterminant& det_up) const {
//   if (singles_alpha_cache_status[i] == CACHED) {
//     return singles_from_alpha.at(i);
//   }

//   std::vector<int> res;

//   const auto& up_elecs = SpinDetUtil::get_occupied_orbitals(det_up);
//   data::SpinDeterminant det_up_copy(det_up);
//   std::unordered_set<int> alpha_singles_set;
//   for (int j = 0; j < n_up; j++) {
//     SpinDetUtil::set_occupation(&det_up_copy, up_elecs[j], false);
//     const auto& alpha_m1 = det_up_copy.SerializeAsString();
//     for (const int alpha_single : unique_ab_m1.at(alpha_m1).first) {
//       if (alpha_singles_set.count(alpha_single) == 0) {
//         res.push_back(alpha_single);
//         alpha_singles_set.insert(alpha_single);
//       }
//     }
//     SpinDetUtil::set_occupation(&det_up_copy, up_elecs[j], true);
//   }
//   std::sort(res.begin(), res.end());

//   return res;
// }

// std::vector<int> ConnectionsSandeepImpl::get_beta_singles(
//     const int i, const data::SpinDeterminant& det_dn) const {
//   if (singles_beta_cache_status[i] == CACHED) {
//     return singles_from_beta.at(i);
//   }

//   std::vector<int> res;

//   const auto& dn_elecs = SpinDetUtil::get_occupied_orbitals(det_dn);
//   data::SpinDeterminant det_dn_copy(det_dn);
//   std::unordered_set<int> beta_singles_set;
//   for (int j = 0; j < n_dn; j++) {
//     SpinDetUtil::set_occupation(&det_dn_copy, dn_elecs[j], false);
//     const auto& beta_m1 = det_dn_copy.SerializeAsString();
//     for (const int beta_single : unique_ab_m1.at(beta_m1).second) {
//       if (beta_singles_set.count(beta_single) == 0) {
//         res.push_back(beta_single);
//         beta_singles_set.insert(beta_single);
//       }
//     }
//     SpinDetUtil::set_occupation(&det_dn_copy, dn_elecs[j], true);
//   }
//   std::sort(res.begin(), res.end());

//   return res;
// }

// void ConnectionsSandeepImpl::sort_by_first(
//     std::vector<int>& vec1, std::vector<int>& vec2) {
//   std::vector<std::pair<int, int>> vec;
//   const int n_vec = vec1.size();
//   for (int i = 0; i < n_vec; i++) {
//     vec.push_back(std::make_pair(vec1[i], vec2[i]));
//   }
//   std::sort(vec.begin(), vec.end(), [&](const auto& a, const auto& b) {
//     return a.first < b.first;
//   });
//   for (int i = 0; i < n_vec; i++) {
//     vec1[i] = vec[i].first;
//     vec2[i] = vec[i].second;
//   }
// }

// std::vector<std::pair<int, double>> ConnectionsSandeepImpl::get_connections(
//     const int i) {
//   if (cache_status[i] == CACHED) {
//     return cached_connections[i];
//   }

//   std::vector<std::pair<int, double>> res;

//   // If already cached in the previous iteration, start from previous
//   results. if (cache_status[i] == CACHE_OUTDATED) {
//     res = std::move(cached_connections[i]);
//   }

//   const auto& det = abstract_system->wf->terms(i).det();

//   // Diagonal.
//   if (cache_status[i] != CACHE_OUTDATED) {
//     // Otherwise, diagonal term must be in the cache already.
//     const double H = abstract_system->hamiltonian(&det, &det);
//     res.push_back(std::make_pair(i, H));
//   }

//   const int start_id = cache_status[i] == CACHE_OUTDATED ? n_dets_prev : i +
//   1;

//   // Two up excitation.
//   const auto& beta = det.dn().SerializeAsString();
//   const int beta_id = unique_betas[beta];
//   const auto& two_ups = beta_major_to_det[beta_id];
//   for (auto it = two_ups.begin(); it != two_ups.end(); it++) {
//     const int det_id = *it;
//     if (det_id < start_id) continue;
//     const auto& det_id_det = abstract_system->wf->terms(det_id).det();
//     const double H = abstract_system->hamiltonian(&det, &det_id_det);
//     if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
//     res.push_back(std::make_pair(det_id, H));
//   }

//   // Two down excitation.
//   const auto& alpha = det.up().SerializeAsString();
//   const int alpha_id = unique_alphas[alpha];
//   const auto& two_dns = alpha_major_to_det[alpha_id];
//   for (auto it = two_dns.begin(); it != two_dns.end(); it++) {
//     const int det_id = *it;
//     if (det_id < start_id) continue;
//     const auto& det_id_det = abstract_system->wf->terms(det_id).det();
//     const double H = abstract_system->hamiltonian(&det, &det_id_det);
//     if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
//     res.push_back(std::make_pair(det_id, H));
//   }

//   // One up one down excitation.
//   if (n_dets - n_dets_prev < 0.2 * n_dets) {
//     std::vector<int> one_up_dets;
//     const auto& single_alphas = get_alpha_singles(alpha_id, det.up());
//     const auto& single_betas = get_beta_singles(beta_id, det.dn());
//     const int thread_id = omp_get_thread_num();
//     for (auto it = single_alphas.begin(); it != single_alphas.end(); it++) {
//       const int single_alpha_id = *it;
//       const auto& alpha_dets = alpha_major_to_det2[single_alpha_id];
//       const auto& start_it =
//           std::lower_bound(alpha_dets.begin(), alpha_dets.end(), start_id);
//       for (auto it_det = start_it; it_det < alpha_dets.end(); it_det++) {
//         const int det_id = *it_det;
//         one_up[thread_id][det_id] = true;
//         one_up_dets.push_back(det_id);
//       }
//     }
//     for (auto it = single_betas.begin(); it != single_betas.end(); it++) {
//       const int single_beta_id = *it;
//       const auto& beta_dets = beta_major_to_det2[single_beta_id];
//       const auto& start_it =
//           std::lower_bound(beta_dets.begin(), beta_dets.end(), start_id);
//       for (auto it_det = start_it; it_det < beta_dets.end(); it_det++) {
//         const int det_id = *it_det;
//         if (one_up[thread_id][det_id]) {
//           const auto& det_id_det = abstract_system->wf->terms(det_id).det();
//           const double H = abstract_system->hamiltonian(&det, &det_id_det);
//           if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
//           res.push_back(std::make_pair(det_id, H));
//         }
//       }
//     }
//     for (const int det_id : one_up_dets) one_up[thread_id][det_id] = false;
//   } else {
//     const auto& one_ups = get_alpha_singles(alpha_id, det.up());
//     const auto& one_dns = get_beta_singles(beta_id, det.dn());
//     for (auto it_up = one_ups.begin(); it_up != one_ups.end(); it_up++) {
//       const int single_alpha_id = *it_up;
//       const auto& connected_beta_ids = alpha_major_to_beta[single_alpha_id];
//       const auto& connected_det_ids = alpha_major_to_det[single_alpha_id];
//       const int n_connected_betas = connected_beta_ids.size();
//       int ptr = 0;
//       for (auto it_dn = one_dns.begin(); it_dn != one_dns.end(); it_dn++) {
//         const int single_beta_id = *it_dn;

//         while (ptr < n_connected_betas &&
//                connected_beta_ids[ptr] < single_beta_id) {
//           ptr++;
//         }
//         if (ptr == n_connected_betas) break;

//         if (connected_beta_ids[ptr] == single_beta_id) {
//           const int det_id = connected_det_ids[ptr];
//           ptr++;
//           if (det_id < start_id) continue;
//           const auto& det_id_det = abstract_system->wf->terms(det_id).det();
//           const double H = abstract_system->hamiltonian(&det, &det_id_det);
//           if (std::abs(H) < std::numeric_limits<double>::epsilon()) continue;
//           res.push_back(std::make_pair(det_id, H));
//         }
//       }
//     }
//   }

//   // Cache if within threshold.
//   if (res.size() <= cache_size) {
//     cached_connections[i] = res;
//     cache_status[i] = CACHED;
//   } else {
//     cached_connections[i].clear();
//     cache_status[i] = CACHE_OVERFLOW;
//   }

//   return res;
// }

// // Connections* Injector::new_connections(
// //     Session* const session, AbstractSystem* const abstract_system) {
// //   return new ConnectionsSandeepImpl(session, abstract_system);
// // }
