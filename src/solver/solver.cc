#include "solver.h"

#include <boost/format.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <unordered_set>
#include <utility>
#include "../injector.h"
#include "../omp_hash_map/src/omp_hash_map.h"
#include "../omp_hash_map/src/reducer.h"
#include "../vector_stats.h"
#include "davidson_util.h"
#include "spin_det_util.h"

#define ENERGY_FORMAT "%.12f"
#define TABLE_FORMAT_ITEM_NAME "%20s"
#define TABLE_FORMAT_F "%17.10f"
#define TABLE_FORMAT_D "%'17d"
#define TABLE_FORMAT_LL "%'17llu"
#define MAX_HASH_LOAD 1.618
#define MAX_DETS_PER_FILE 20000000

class UncertainResult {
 public:
  double value = 0.0;
  double uncertainty = 0.0;
};

class SolverImpl : public Solver {
 public:
  ~SolverImpl();

  SolverImpl(
      Session* const session,
      Connections* const connections,
      AbstractSystem* const abstrct_system);

  void setup_hf() override;

  void variation(const double eps_var) override;

  void save_variation_result(const std::string& filename) override;

  bool load_variation_result(const std::string& filename) override;

  std::vector<double> apply_hamiltonian(
      const std::vector<double>& vec, const bool first_iteration) override;

  void perturbation(
      const int n_orbs_var,
      const double eps_var,
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts) override;

 private:
  int n_up = 0;

  int n_dn = 0;

  double energy_hf = 0.0;

  double energy_var = 0.0;

  std::unordered_set<std::string> var_dets_set;

  bool verbose = false;

  void print_var_result() const;

  std::vector<data::Determinant> tmp_dets;

  std::vector<double> get_energy_pts_pre_dtm(
      const std::vector<int>& n_orbs_pts);

  std::vector<UncertainResult> get_energy_pts_dtm(
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& energy_pts_pre_dtm);

  std::vector<std::vector<UncertainResult>> get_energy_pts_stc(
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts,
      const std::vector<UncertainResult>& energy_pts_dtm);

  double get_weight(const int i) const;
};

SolverImpl::SolverImpl(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system)
    : Solver(session, connections, abstract_system) {
  Parallel* const parallel = session->get_parallel();
  Config* const config = session->get_config();
  verbose = parallel->is_master();
  tmp_dets.resize(parallel->get_n_threads());
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
}

SolverImpl::~SolverImpl() {}

void SolverImpl::setup_hf() {
  // Create or clear wavefunction.
  abstract_system->dets.clear();
  abstract_system->coefs.clear();

  // Add a single term with coef 1.0 and no diffs.
  auto& det_hf = tmp_dets[0];
  det_hf.mutable_up()->set_n_hf_elecs(n_up);
  det_hf.mutable_dn()->set_n_hf_elecs(n_dn);
  abstract_system->dets.push_back(det_hf.SerializeAsString());
  abstract_system->coefs.push_back(1.0);

  // Update HF and variational energy.
  energy_hf = energy_var = abstract_system->hamiltonian(&det_hf, &det_hf);

  // Clear connections.
  connections->clear();

  if (verbose) printf("HF det setup.\n");
};

void SolverImpl::variation(const double eps_var) {
  Timer* const timer = session->get_timer();

  // Contruct variational determinant hash set.
  const int n_start_dets = abstract_system->dets.size();
  var_dets_set.clear();
  std::vector<double> prev_coefs;
  var_dets_set.reserve(n_start_dets);
  prev_coefs.resize(n_start_dets, 0.0);
  for (int i = 0; i < n_start_dets; i++) {
    var_dets_set.insert(abstract_system->dets[i]);
  }

  // Variation iterations.
  bool converged = false;
  int iteration = 0;
  int n_new_dets = 0;

  // Add to the var wf if not exists and set its coef to zero.
  const auto& connected_det_handler =
      [&](const data::Determinant* const connected_det) {
        std::string det_code = connected_det->SerializeAsString();
        if (var_dets_set.count(det_code) == 0) {
          var_dets_set.insert(det_code);
          abstract_system->dets.push_back(det_code);
          abstract_system->coefs.push_back(0.0);
          n_new_dets++;
        }
      };

  // Apply hamiltonian to a vec (with partially cached hamiltonian).
  const auto& apply_hamiltonian_func = std::bind(
      &SolverImpl::apply_hamiltonian,
      this,
      std::placeholders::_1,
      std::placeholders::_2);

  while (!converged) {
    timer->start("loop " + std::to_string(iteration));

    // Find new dets.
    n_new_dets = 0;
    const int n_old_dets = abstract_system->dets.size();
    for (int i = 0; i < n_old_dets; i++) {
      const double coef = abstract_system->coefs[i];
      if (i == 0) printf("coef: %f\n", coef);
      auto& det = tmp_dets[0];
      det.ParseFromString(abstract_system->dets[i]);
      if (std::abs(coef) <= std::abs(prev_coefs[i])) continue;
      abstract_system->find_connected_dets(
          &det, eps_var / std::abs(coef), connected_det_handler);
    }
    const int n_total_dets = abstract_system->dets.size();
    if (verbose) {
      printf("New / total dets: %'d / %'d\n", n_new_dets, n_total_dets);
    }
    timer->checkpoint("found new dets");

    // Update prev coefs.
    prev_coefs.resize(n_total_dets);
    for (int i = 0; i < n_total_dets; i++) {
      prev_coefs[i] = abstract_system->coefs[i];
    }

    // Diagonalize.
    connections->update();
    timer->checkpoint("updated connections");
    std::vector<double> diagonal(n_total_dets, 0.0);
#pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < n_total_dets; i++) {
      const int thread_id = omp_get_thread_num();
      auto& tmp_det = tmp_dets[thread_id];
      tmp_det.ParseFromString(abstract_system->dets[i]);
      diagonal[i] = abstract_system->hamiltonian(&tmp_det, &tmp_det);
    }
    const auto& diagonalization_result = DavidsonUtil::diagonalize(
        prev_coefs, diagonal, apply_hamiltonian_func, 5, verbose);
    const double energy_var_new = diagonalization_result.first;
    const auto& new_coefs = diagonalization_result.second;
    for (int i = 0; i < n_total_dets; i++) {
      abstract_system->coefs[i] = new_coefs[i];
    }

    // Determine convergence.
    if (std::abs(energy_var_new - energy_var) < 1.0e-7) {
      converged = true;
    }

    energy_var = energy_var_new;
    iteration++;

    if (verbose) print_var_result();

    timer->end();  // iteration.
  }
};

void SolverImpl::save_variation_result(const std::string& filename) {
  // Write summary.
  std::fstream var_file(
      filename, std::ios::out | std::ios::trunc | std::ios::binary);
  data::VariationResult res;
  data::Wavefunction trunk_wf;
  res.set_energy_hf(energy_hf);
  res.set_energy_var(energy_var);
  const int n_dets = abstract_system->dets.size();
  res.set_n_dets(n_dets);
  res.SerializeToOstream(&var_file);
  var_file.close();
  if (verbose) {
    printf("Saved summary to: %s\n", filename.c_str());
  }

  // Write wavefunction trunks. (Due to the 2GB serialization limit of Protobuf)
  int trunk_n_dets = 0;
  int trunk_id = 0;
  for (int i = 0; i < n_dets; i++) {
    data::Term* const new_term = trunk_wf.add_terms();
    new_term->set_coef(abstract_system->coefs[i]);
    new_term->mutable_det()->ParseFromString(abstract_system->dets[i]);
    trunk_n_dets++;
    if (trunk_n_dets >= MAX_DETS_PER_FILE || i == n_dets - 1) {
      const auto& wf_filename = filename + ".part" + std::to_string(trunk_id);
      std::fstream wf_file(
          wf_filename, std::ios::out | std::ios::trunc | std::ios::binary);
      trunk_wf.SerializeToOstream(&wf_file);
      wf_file.close();
      if (verbose) {
        printf("Saved %'d dets to: %s\n", trunk_n_dets, wf_filename.c_str());
      }
      trunk_wf.Clear();
      trunk_n_dets = 0;
      trunk_id++;
    }
  }
}

bool SolverImpl::load_variation_result(const std::string& filename) {
  // Load summary.
  std::fstream var_file(filename, std::ios::in | std::ios::binary);
  data::VariationResult res;
  if (!res.ParseFromIstream(&var_file)) {
    return false;
  }
  var_file.close();
  energy_hf = res.energy_hf();
  energy_var = res.energy_var();
  const int n_dets_total = res.n_dets();
  if (verbose) {
    printf("Loaded summary from: %s\n", filename.c_str());
  }

  // Load wavefunctions from trunks.
  abstract_system->dets.clear();
  abstract_system->coefs.clear();
  int n_dets_read = 0;
  int trunk_id = 0;
  data::Wavefunction trunk_wf;
  while (n_dets_read < n_dets_total) {
    const auto& wf_filename = filename + ".part" + std::to_string(trunk_id);
    trunk_id++;
    std::fstream var_file(wf_filename, std::ios::in | std::ios::binary);
    if (!trunk_wf.ParseFromIstream(&var_file)) {
      throw std::runtime_error("variational results corrupted");
    }
    const int trunk_n_dets = trunk_wf.terms_size();
    for (int i = 0; i < trunk_n_dets; i++) {
      const auto& term = trunk_wf.terms(i);
      abstract_system->coefs.push_back(term.coef());
      abstract_system->dets.push_back(term.det().SerializeAsString());
      n_dets_read++;
    }
    if (verbose) {
      printf("Loaded %'d dets from: %s\n", trunk_n_dets, wf_filename.c_str());
    }
  }

  if (verbose) {
    print_var_result();
  }

  return true;
}

void SolverImpl::print_var_result() const {
  printf("Number of dets: %'zu\n", abstract_system->dets.size());
  printf("Variation energy: " ENERGY_FORMAT " Ha\n", energy_var);
  const double energy_corr = energy_var - energy_hf;
  printf("Correlation energy (variation): " ENERGY_FORMAT " Ha\n", energy_corr);
}

std::vector<double> SolverImpl::apply_hamiltonian(
    const std::vector<double>& vec, const bool first_iteration) {
  const int n_dets = vec.size();
  Parallel* const parallel = session->get_parallel();
  const int proc_id = parallel->get_proc_id();
  const int n_procs = parallel->get_n_procs();
  const int n_threads = parallel->get_n_threads();
  std::vector<std::vector<double>> res(n_threads);
  for (int i = 0; i < n_threads; i++) res[i].resize(n_dets, 0.0);
  std::vector<unsigned long long> n_nonzero_elems(n_threads, 0);

#pragma omp parallel for
  for (int i = proc_id; i < n_dets; i += n_procs) {
    const int thread_id = omp_get_thread_num();
    const auto& conns = connections->get_connections(i);
    for (const auto conn : conns) {
      const int j = conn.first;
      const double H_ij = conn.second;
      res[thread_id][i] += H_ij * vec[j];
      if (i != j) {
        res[thread_id][j] += H_ij * vec[i];
        n_nonzero_elems[thread_id] += 2;
      } else {
        n_nonzero_elems[thread_id]++;
      }
    }
  }

  for (int i = 1; i < n_threads; i++) {
    n_nonzero_elems[0] += n_nonzero_elems[i];
    for (int j = 0; j < n_dets; j++) {
      res[0][j] += res[i][j];
    }
  }
  parallel->reduce_to_sum(res[0]);
  parallel->reduce_to_sum(n_nonzero_elems[0]);
  if (verbose && first_iteration) {
    printf("Number of non-zero elements: %'llu\n", n_nonzero_elems[0]);
  }

  session->get_timer()->checkpoint("hamiltonian applied");

  return res[0];
};

void SolverImpl::perturbation(
    const int n_orbs_var,
    const double eps_var,
    const std::vector<int>& n_orbs_pts,
    const std::vector<double>& eps_pts) {
  // Check if the results already exists.
  const auto result_filename =
      str(boost::format("pt_%d_%#.4g.csv") % n_orbs_var % eps_var);
  if (std::ifstream(result_filename)) {
    if (verbose) printf("PT results found in: %s\n", result_filename.c_str());
    return;
  }

  // Clean variation variables.
  connections->clear();

  // Construct var dets set.
  var_dets_set.clear();
  var_dets_set.reserve(abstract_system->dets.size());
  for (const auto& det_string : abstract_system->dets) {
    var_dets_set.insert(det_string);
  }

  // Get Deterministic PT correction.
  const auto& energy_pts_pre_dtm = get_energy_pts_pre_dtm(n_orbs_pts);
  const auto& energy_pts_dtm =
      get_energy_pts_dtm(n_orbs_pts, energy_pts_pre_dtm);

  // Check eps_pts validity.
  const double eps_dtm_pt = config->get_double("eps_dtm_pt");
  const int n_eps_pts = eps_pts.size();
  for (int i = 1; i < n_eps_pts; i++) {
    assert(eps_pts[i] <= eps_dtm_pt);
    assert(eps_pts[i] < eps_pts[i - 1]);
  }

  // Get Stochastic PT correction.
  const auto& energy_pts_stc =
      get_energy_pts_stc(n_orbs_pts, eps_pts, energy_pts_dtm);

  // Record results.
  std::ofstream result_file(result_filename);
  result_file << "n_orbs_var,eps_var,n_orbs_pt,eps_pt,energy_corr,uncert,"
              << "energy_hf,energy_var,energy_pt" << std::endl;

  const int n_n_orbs_pts = n_orbs_pts.size();
  for (int i = 0; i < n_eps_pts; i++) {
    for (int j = 0; j < n_n_orbs_pts; j++) {
      const double eps_pt = eps_pts[i];
      const int n_orbs_pt = n_orbs_pts[j];
      const double energy_pt_dtm = energy_pts_dtm[j].value;
      const double energy_pt_stc = energy_pts_stc[i][j].value;
      const double energy_pt = energy_pt_dtm + energy_pt_stc;
      const double uncert = std::sqrt(
          std::pow(energy_pts_stc[i][j].uncertainty, 2) +
          std::pow(energy_pts_dtm[j].uncertainty, 2));
      const double energy_corr = energy_var + energy_pt - energy_hf;
      result_file << str(boost::format("%d, %#.4g, %d, %#.4g, %#.15g, %#.15g, "
                                       "%#.15g, %#.15g, %#.15g") %
                         n_orbs_var % eps_var % n_orbs_pt % eps_pt %
                         energy_corr % uncert % energy_hf % energy_var %
                         energy_pt)
                  << std::endl;
    }
  }

  if (verbose) {
    printf("Results saved to: %s\n", result_filename.c_str());
  }

  result_file.close();
}

std::vector<double> SolverImpl::get_energy_pts_pre_dtm(
    const std::vector<int>& n_orbs_pts) {
  // Cache commonly used variables.
  const int n_n_orbs_pts = n_orbs_pts.size();
  const int n_var_dets = var_dets_set.size();
  const double eps_pre_dtm_pt = config->get_double("eps_pre_dtm_pt");
  std::hash<std::string> string_hasher;

  // Construct partial sums store and pt results store.
  omp_hash_map<std::string, double> partial_sums_pre;
  partial_sums_pre.set_max_load_factor(MAX_HASH_LOAD);
  std::vector<double> energy_pts_pre_dtm(n_n_orbs_pts, 0.0);
  std::vector<unsigned long long> n_pt_dets_pre_dtm(n_n_orbs_pts, 0);

  timer->start("pre_dtm");
  if (verbose) printf(">>> eps_pre_dtm_pt %#.4g\n", eps_pre_dtm_pt);
  timer->start("search");
  double target_progress = 0.25;
#pragma omp parallel for schedule(dynamic, 5)
  for (int i = 0; i < n_var_dets; i++) {
    auto& var_det = tmp_dets[omp_get_thread_num()];
    var_det.ParseFromString(abstract_system->dets[i]);
    const double coef = abstract_system->coefs[i];

    const auto& pt_det_handler = [&](const auto& det_a) {
      const auto& det_a_code = det_a->SerializeAsString();
      const size_t det_a_hash = string_hasher(det_a_code);
      if (var_dets_set.count(det_a_code) == 1) return;
      if (det_a_hash % n_procs != proc_id) return;
      const double H_ai = abstract_system->hamiltonian(&var_det, det_a);
      const double partial_sum_term = H_ai * coef;
      partial_sums_pre.set(
          det_a_code, [&](double& value) { value += partial_sum_term; }, 0.0);
    };

    abstract_system->find_connected_dets(
        &var_det, eps_pre_dtm_pt / std::abs(coef), pt_det_handler);

    // Report progress. (every time progress reaches 2^n %).
    const double current_progress = i * 100.0 / n_var_dets;
    if (target_progress <= current_progress) {
      if (omp_get_thread_num() == 0) {
        if (verbose) {
          printf("Master node PT dets: %'zu\n", partial_sums_pre.get_n_keys());
        }
        std::string event =
            str(boost::format("Progress: %.2f %%") % target_progress);
        timer->checkpoint(event);
        target_progress *= 2.0;
      }
    }
  }
  while (target_progress <= 100) {
    std::string event =
        str(boost::format("Progress: %.2f %%") % target_progress);
    timer->checkpoint(event);
    target_progress *= 2.0;
  }
  timer->end();  // Search.

  timer->start("accumulate");
  partial_sums_pre.apply([&](const std::string& det_code, const double value) {
    const int thread_id = omp_get_thread_num();
    auto& det = tmp_dets[thread_id];
    det.ParseFromString(det_code);
    const int n_orbs_used = SpinDetUtil::get_n_orbs_used(det);
    const double H_aa = abstract_system->hamiltonian(&det, &det);
    const double contribution = value * value / (energy_var - H_aa);
    for (int i = 0; i < n_n_orbs_pts; i++) {
      if (n_orbs_used > n_orbs_pts[i]) continue;
#pragma omp atomic
      n_pt_dets_pre_dtm[i]++;
#pragma omp atomic
      energy_pts_pre_dtm[i] += contribution;
    }
  });
  partial_sums_pre.clear();
  parallel->reduce_to_sum(n_pt_dets_pre_dtm);
  parallel->reduce_to_sum(energy_pts_pre_dtm);

  if (verbose) {
    printf(TABLE_FORMAT_ITEM_NAME, "# orbitals PT:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      printf(TABLE_FORMAT_D, n_orbs_pts[i]);
    }
    printf("\n");
    printf(TABLE_FORMAT_ITEM_NAME, "Pre DTM energy PT:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      printf(TABLE_FORMAT_F, energy_pts_pre_dtm[i]);
    }
    printf("\n");
    printf(TABLE_FORMAT_ITEM_NAME, "EST. CORR. energy:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      const double energy_corr = energy_pts_pre_dtm[i] + energy_var - energy_hf;
      printf(TABLE_FORMAT_F, energy_corr);
    }
    printf("\n");
    printf(TABLE_FORMAT_ITEM_NAME, "# DTM PT dets:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      printf(TABLE_FORMAT_LL, n_pt_dets_pre_dtm[i]);
    }
    printf("\n");
  }
  timer->end();  // Accumulate.
  timer->end();  // pre dtm.

  return energy_pts_pre_dtm;
}

std::vector<UncertainResult> SolverImpl::get_energy_pts_dtm(
    const std::vector<int>& n_orbs_pts,
    const std::vector<double>& energy_pts_pre_dtm) {
  // Cache commonly used variables.
  const int n_n_orbs_pts = n_orbs_pts.size();
  const int n_var_dets = var_dets_set.size();
  const size_t n_pt_batches_dtm = config->get_int("n_batches_dtm_pt");
  const double eps_pre_dtm_pt = config->get_double("eps_pre_dtm_pt");
  const double eps_dtm_pt = config->get_double("eps_dtm_pt");
  std::hash<std::string> string_hasher;

  // Construct partial sums store and pt results store.
  omp_hash_map<std::string, double> partial_sums_pre;
  omp_hash_map<std::string, double> partial_sums;
  partial_sums.set_max_load_factor(MAX_HASH_LOAD);
  partial_sums_pre.set_max_load_factor(MAX_HASH_LOAD);
  std::vector<UncertainResult> energy_pts_dtm(n_n_orbs_pts);
  std::vector<std::vector<double>> energy_pts_dtm_batches(n_n_orbs_pts);
  std::vector<unsigned long long> n_pt_dets_dtm(n_n_orbs_pts, 0);
  const double target_error = config->get_double("target_error");

  timer->start("dtm");
  if (verbose) printf(">>> eps_dtm_pt %#.4g\n", eps_dtm_pt);
  // Process batch by batch to achieve a larger run in a constrained mem env.
  for (size_t b = 0; b < n_pt_batches_dtm; b++) {
    timer->start(str(boost::format("%d/%d") % (b + 1) % n_pt_batches_dtm));

    double target_progress = 0.25;
    timer->start("search");
#pragma omp parallel for schedule(dynamic, 5)
    for (int i = 0; i < n_var_dets; i++) {
      auto& var_det = tmp_dets[omp_get_thread_num()];
      var_det.ParseFromString(abstract_system->dets[i]);
      const double coef = abstract_system->coefs[i];

      // Process only if:
      // 1. is not a var det, and
      // 2. belongs to the current proc_id, and
      // 3. belongs to the current batch.
      const auto& pt_det_handler = [&](const auto& det_a) {
        const auto& det_a_code = det_a->SerializeAsString();
        const size_t det_a_hash = string_hasher(det_a_code);
        if (var_dets_set.count(det_a_code) == 1) return;
        if (det_a_hash % n_procs != proc_id) return;
        if ((det_a_hash / n_procs) % n_pt_batches_dtm != b) return;
        const double H_ai = abstract_system->hamiltonian(&var_det, det_a);
        const double partial_sum_term = H_ai * coef;
        if (std::abs(partial_sum_term) >= eps_pre_dtm_pt) {
          partial_sums_pre.set(
              det_a_code,
              [&](double& value) { value += partial_sum_term; },
              0.0);
        }
        partial_sums.set(
            det_a_code, [&](double& value) { value += partial_sum_term; }, 0.0);
      };

      abstract_system->find_connected_dets(
          &var_det, eps_dtm_pt / std::abs(coef), pt_det_handler);

      // Report progress. (every time progress reaches 2^n %).
      const double current_progress = i * 100.0 / n_var_dets;
      if (target_progress <= current_progress) {
        if (omp_get_thread_num() == 0) {
          if (verbose) {
            printf("Master node PT dets: %'zu\n", partial_sums.get_n_keys());
          }
          std::string event =
              str(boost::format("Progress: %.2f %%") % target_progress);
          timer->checkpoint(event);
          target_progress *= 2.0;
        }
      }
    }
    while (target_progress <= 100) {
      std::string event =
          str(boost::format("Progress: %.2f %%") % target_progress);
      timer->checkpoint(event);
      target_progress *= 2.0;
    }
    timer->end();  // Search.

    timer->start("accumulate");

    // For reduction after each batch.
    std::vector<double> energy_pts_dtm_batch(n_n_orbs_pts, 0.0);

    // Apply to each key, value pair:
    // 1. Calculate the contribution.
    // 2. For each n_orbs_pt interested:
    //      If the n_orbs used by the key is less than n_orbs_pt:
    //        Increase n_pt_dets by one and add contribution to energy_pt.
    partial_sums.apply([&](const std::string& det_code, const double value) {
      const int thread_id = omp_get_thread_num();
      auto& det = tmp_dets[thread_id];
      det.ParseFromString(det_code);
      const int n_orbs_used = SpinDetUtil::get_n_orbs_used(det);
      const double H_aa = abstract_system->hamiltonian(&det, &det);
      const double value_pre =
          partial_sums_pre.get_copy_or_default(det_code, 0.0);
      partial_sums_pre.unset(det_code);
      const double contribution =
          (value * value - value_pre * value_pre) / (energy_var - H_aa);
      for (int i = 0; i < n_n_orbs_pts; i++) {
        if (n_orbs_used > n_orbs_pts[i]) continue;
#pragma omp atomic
        n_pt_dets_dtm[i]++;
#pragma omp atomic
        energy_pts_dtm_batch[i] += contribution;
      }
    });

    // Aggregate the results from each proc.
    parallel->reduce_to_sum(energy_pts_dtm_batch);
    double max_uncert = 0.0;
    for (int i = 0; i < n_n_orbs_pts; i++) {
      energy_pts_dtm_batches[i].push_back(energy_pts_dtm_batch[i]);
    }
    for (int i = 0; i < n_n_orbs_pts; i++) {
      energy_pts_dtm[i].value =
          get_avg(energy_pts_dtm_batches[i]) * n_pt_batches_dtm +
          energy_pts_pre_dtm[i];
      if (b == n_pt_batches_dtm - 1) {
        energy_pts_dtm[i].uncertainty = 0.0;  // All batches finished.
      } else {
        energy_pts_dtm[i].uncertainty =
            get_stdev(energy_pts_dtm_batches[i]) * n_pt_batches_dtm / b;
      }
      max_uncert = std::max(max_uncert, energy_pts_dtm[i].uncertainty);
    }

    // Print batch result and estimated total correction.
    if (verbose) {
      printf(TABLE_FORMAT_ITEM_NAME, "# orbitals PT:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf(TABLE_FORMAT_D, n_orbs_pts[i]);
      }
      printf("\n");
      printf(TABLE_FORMAT_ITEM_NAME, "Batch DTM energy PT:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf(TABLE_FORMAT_F, energy_pts_dtm_batch[i]);
      }
      printf("\n");
      printf(TABLE_FORMAT_ITEM_NAME, "EST. DTM energy PT:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf(TABLE_FORMAT_F, energy_pts_dtm[i].value - energy_pts_pre_dtm[i]);
      }
      printf("\n");
      printf(TABLE_FORMAT_ITEM_NAME, "(including pre dtm):");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf(TABLE_FORMAT_F, energy_pts_dtm[i].value);
      }
      printf("\n");
      printf(TABLE_FORMAT_ITEM_NAME, "Uncertainty:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf(TABLE_FORMAT_F, energy_pts_dtm[i].uncertainty);
      }
      printf("\n");
      printf(TABLE_FORMAT_ITEM_NAME, "EST. CORR. energy:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        const double energy_corr =
            energy_pts_dtm[i].value + energy_var - energy_hf;
        printf(TABLE_FORMAT_F, energy_corr);
      }
      printf("\n");
    }

    timer->end();  // Accumulation.

    partial_sums.clear();
    partial_sums_pre.clear();
    timer->end();  // Batch.

    if (b > 0 && b < n_pt_batches_dtm - 2 && max_uncert < 0.2 * target_error) {
      if (verbose) {
        printf(
            "\n>>> Skip remaining batches since uncertainty is significantly "
            "smaller than target error.\n");
      }
      break;
    }
  }

  parallel->reduce_to_sum(n_pt_dets_dtm);

  // Print total PT dets.
  if (verbose) {
    printf(TABLE_FORMAT_ITEM_NAME, "# DTM PT dets:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      printf(TABLE_FORMAT_LL, n_pt_dets_dtm[i]);
    }
    printf("\n");
  }

  timer->end();

  return energy_pts_dtm;
}

std::vector<std::vector<UncertainResult>> SolverImpl::get_energy_pts_stc(
    const std::vector<int>& n_orbs_pts,
    const std::vector<double>& eps_pts,
    const std::vector<UncertainResult>& energy_pts_dtm) {
  const int n_n_orbs_pts = n_orbs_pts.size();
  const int n_eps_pts = eps_pts.size();

  // Construct results store.
  std::vector<std::vector<UncertainResult>> energy_pts_stc(n_eps_pts);
  std::vector<std::vector<std::vector<double>>> energy_pts_loops(n_eps_pts);
  std::vector<omp_hash_map<std::string, double>> partial_sums(n_eps_pts);
  omp_hash_map<std::string, double> partial_sums_dtm;
  for (int i = 0; i < n_eps_pts; i++) {
    energy_pts_stc[i].resize(n_n_orbs_pts);
    energy_pts_loops[i].resize(n_n_orbs_pts);
    partial_sums[i].set_max_load_factor(MAX_HASH_LOAD);
  }
  partial_sums_dtm.set_max_load_factor(MAX_HASH_LOAD);

  // Construct probabilities from normalized weights.
  const int n_var_dets = var_dets_set.size();
  double sum_weights = 0.0;
  for (int i = 0; i < n_var_dets; i++) {
    sum_weights += get_weight(i);
  }
  std::vector<double> probs(n_var_dets);
  std::vector<double> cum_probs(n_var_dets);  // For sampling.
  for (int i = 0; i < n_var_dets; i++) {
    probs[i] = get_weight(i) / sum_weights;
    if (i == 0) {
      cum_probs[i] = probs[i];
    } else {
      cum_probs[i] = probs[i] + cum_probs[i - 1];
    }
  }

  // Stochastic iterations.
  int iteration = 0;
  const int max_n_iterations = config->get_int("max_n_iterations_stc_pt");
  std::unordered_map<int, int> stc_pt_sample_dets;
  std::vector<int> stc_pt_sample_dets_list;
  const int n_samples = config->get_int("n_samples_stc_pt");
  srand(time(NULL));
  const double eps_dtm_pt = config->get_double("eps_dtm_pt");
  const double eps_pts_min = eps_pts.back();
  const size_t n_batches = config->get_int("n_batches_stc_pt");
  std::hash<std::string> string_hasher;
  const double target_error = config->get_double("target_error");

  timer->start("stc");
  while (iteration < max_n_iterations) {
    timer->start(str(boost::format("loop %d") % iteration));

    // Cleanup.
    for (int i = 0; i < n_eps_pts; i++) {
      partial_sums[i].clear();
    }
    partial_sums_dtm.clear();
    stc_pt_sample_dets.clear();
    stc_pt_sample_dets_list.clear();

    std::vector<std::vector<double>> energy_pts_loop(n_eps_pts);
    std::vector<std::vector<unsigned long long>> n_pt_dets_loop(n_eps_pts);
    for (int i = 0; i < n_eps_pts; i++) {
      energy_pts_loop[i].resize(n_n_orbs_pts, 0.0);
      n_pt_dets_loop[i].resize(n_n_orbs_pts, 0);
    }

    // Generate random samples.
    for (int i = 0; i < n_samples; i++) {
      const double rand_01 = ((double)rand() / (RAND_MAX));
      const int sample_det_id =
          std::lower_bound(cum_probs.begin(), cum_probs.end(), rand_01) -
          cum_probs.begin();
      if (stc_pt_sample_dets.count(sample_det_id) == 0) {
        stc_pt_sample_dets_list.push_back(sample_det_id);
      }
      stc_pt_sample_dets[sample_det_id]++;
    }

    // Use one batch to estimate contribution from all batches.
    const size_t selected_batch = rand() % n_batches;

    // Search PT dets from selected sample var dets.
    const int n_stc_pt_sample_dets = stc_pt_sample_dets_list.size();

#pragma omp parallel for schedule(dynamic, 2)
    for (int s = 0; s < n_stc_pt_sample_dets; s++) {
      const int var_det_id = stc_pt_sample_dets_list[s];
      const double prob = probs[var_det_id];
      const double cnt = static_cast<double>(stc_pt_sample_dets[var_det_id]);
      auto& var_det = tmp_dets[omp_get_thread_num()];
      var_det.ParseFromString(abstract_system->dets[var_det_id]);
      const double coef = abstract_system->coefs[var_det_id];

      const auto& pt_det_handler = [&](const auto& det_a) {
        const auto& det_a_code = det_a->SerializeAsString();
        const size_t det_a_hash = string_hasher(det_a_code);
        if (var_dets_set.count(det_a_code) == 1) return;
        if (det_a_hash % n_procs != proc_id) return;
        if ((det_a_hash / n_procs) % n_batches != selected_batch) {
          return;  // Only use one batch.
        }
        const double H_ai = abstract_system->hamiltonian(&var_det, det_a);
        const double H_aa = abstract_system->hamiltonian(det_a, det_a);
        const double factor =
            n_batches / ((energy_var - H_aa) * n_samples * (n_samples - 1));
        const double partial_sum_term = H_ai * coef;
        const double contrib_1 = cnt * partial_sum_term / prob * sqrt(-factor);

        // Add to dtm partial sums.
        if (std::abs(partial_sum_term) >= eps_dtm_pt) {
          partial_sums_dtm.set(
              det_a_code, [&](double& value) { value += contrib_1; }, 0.0);
        }

        // Add to partial sums.
        for (int i = 0; i < n_eps_pts; i++) {
          if (std::abs(partial_sum_term) >= eps_pts[i]) {
            partial_sums[i].set(
                det_a_code, [&](double& value) { value += contrib_1; }, 0.0);
            break;
          }
        }

        // Calculate 2nd order contribution if not belong to dtm pt .
        if (std::abs(partial_sum_term) >= eps_dtm_pt) return;

        const double contrib_2 =
            (cnt * (n_samples - 1) / prob - (cnt * cnt) / (prob * prob)) *
            partial_sum_term * partial_sum_term * factor;
        const int n_orbs_used = SpinDetUtil::get_n_orbs_used(*det_a);
        for (int i = 0; i < n_eps_pts; i++) {
          if (std::abs(partial_sum_term) < eps_pts[i]) continue;
          for (int j = 0; j < n_n_orbs_pts; j++) {
            if (n_orbs_pts[j] < n_orbs_used) continue;
#pragma omp atomic
            energy_pts_loop[i][j] += contrib_2;
          }
        }
      };

      abstract_system->find_connected_dets(
          &var_det, eps_pts_min / std::abs(coef), pt_det_handler);
    }

    // Accumulate and get PT correction from batch.
    for (int i = 0; i < n_eps_pts; i++) {
      partial_sums[i].apply(
          [&](const std::string& det_code, const double value) {
            const int thread_id = omp_get_thread_num();
            auto& det = tmp_dets[thread_id];
            det.ParseFromString(det_code);

            double cum_value = value;
            const double value_dtm =
                partial_sums_dtm.get_copy_or_default(det_code, 0.0);
            partial_sums_dtm.unset(det_code);

            const int n_orbs_used = SpinDetUtil::get_n_orbs_used(det);
            for (int j = i; j < n_eps_pts; j++) {
              if (j > i) {
                const double value_j =
                    partial_sums[j].get_copy_or_default(det_code, 0.0);
                partial_sums[j].unset(det_code);
                cum_value += value_j;
              }

              const double contrib_1 =
                  -cum_value * cum_value + value_dtm * value_dtm;

              for (int k = 0; k < n_n_orbs_pts; k++) {
                if (n_orbs_pts[k] < n_orbs_used) continue;
#pragma omp atomic
                energy_pts_loop[j][k] += contrib_1;
#pragma omp atomic
                n_pt_dets_loop[j][k]++;
              }
            }
          });
    }

    // Append loop results store.
    for (int i = 0; i < n_eps_pts; i++) {
      parallel->reduce_to_sum(energy_pts_loop[i]);
      parallel->reduce_to_sum(n_pt_dets_loop[i]);
      for (int j = 0; j < n_n_orbs_pts; j++) {
        energy_pts_loops[i][j].push_back(energy_pts_loop[i][j]);
      }
    }

    // Calculate statistics.
    double max_uncert = 0.0;
    for (int i = 0; i < n_eps_pts; i++) {
      for (int j = 0; j < n_n_orbs_pts; j++) {
        const double avg = get_avg(energy_pts_loops[i][j]);
        const double stdev = get_stdev(energy_pts_loops[i][j]);
        energy_pts_stc[i][j].value = avg;
        energy_pts_stc[i][j].uncertainty = stdev / sqrt(iteration);
        max_uncert = std::max(energy_pts_stc[i][j].uncertainty, max_uncert);
      }
    }

    // Report results.
    if (verbose) {
      printf(TABLE_FORMAT_ITEM_NAME, "# orbitals PT:");
      for (int j = 0; j < n_n_orbs_pts; j++) {
        printf(TABLE_FORMAT_D, n_orbs_pts[j]);
      }
      printf("\n");
      for (int i = 0; i < n_eps_pts; i++) {
        printf(">>> eps_stc_pt %#.4g\n", eps_pts[i]);
        printf(TABLE_FORMAT_ITEM_NAME, "Loop # STC PT dets:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          printf(TABLE_FORMAT_LL, n_pt_dets_loop[i][j]);
        }
        printf("\n");
        printf(TABLE_FORMAT_ITEM_NAME, "Loop STC energy PT:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          printf(TABLE_FORMAT_F, energy_pts_loops[i][j][iteration]);
        }
        printf("\n");
        printf(TABLE_FORMAT_ITEM_NAME, "EST. STC energy PT:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          printf(TABLE_FORMAT_F, energy_pts_stc[i][j].value);
        }
        printf("\n");
        printf(TABLE_FORMAT_ITEM_NAME, "Uncertainty:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          printf(TABLE_FORMAT_F, energy_pts_stc[i][j].uncertainty);
        }
        printf("\n");
        printf(TABLE_FORMAT_ITEM_NAME, "EST. energy PT:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          printf(
              TABLE_FORMAT_F,
              energy_pts_stc[i][j].value + energy_pts_dtm[j].value);
        }
        printf("\n");
        printf(TABLE_FORMAT_ITEM_NAME, "EST. CORR. energy:");
        for (int j = 0; j < n_n_orbs_pts; j++) {
          const double energy_corr = energy_pts_stc[i][j].value +
                                     energy_pts_dtm[j].value + energy_var -
                                     energy_hf;
          printf(TABLE_FORMAT_F, energy_corr);
        }
        printf("\n");
      }
    }

    timer->end();

    if (iteration >= 10 && max_uncert <= target_error) break;

    iteration++;
  }

  timer->end();

  return energy_pts_stc;
}

double SolverImpl::get_weight(const int i) const {
  const double coef = abstract_system->coefs[i];
  const double abs_coef = std::abs(coef);
  return abs_coef;
}

Solver* Injector::new_solver(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system) {
  return new SolverImpl(session, connections, abstract_system);
}
