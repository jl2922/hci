#include "solver.h"

#include <boost/format.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <unordered_set>
#include <utility>
#include "../injector.h"
#include "../omp_hash_map/src/omp_hash_map.h"
#include "../omp_hash_map/src/reducer.h"
#include "connections.h"
#include "davidson_util.h"
#include "spin_det_util.h"

#define ENERGY_FORMAT "%.12f"

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
      const std::vector<double>& vec) override;

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

  bool verbose = false;

  void print_var_result() const;

  std::vector<double> get_energy_pts_dtm(
      const std::vector<int>& n_orbs_pts) const;

  std::vector<std::vector<double>> get_energy_pts_stc(
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts,
      const std::vector<double>& energy_pts_dtm) const;
};

SolverImpl::SolverImpl(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system)
    : Solver(session, connections, abstract_system) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
}

SolverImpl::~SolverImpl() {}

void SolverImpl::setup_hf() {
  // Create or clear wavefunction.
  if (!abstract_system->wf) {
    abstract_system->wf.reset(new data::Wavefunction());
  }
  data::Wavefunction* const wf = abstract_system->wf.get();
  wf->Clear();

  // Add a single term with coef 1.0 and no diffs.
  data::Term* term = wf->add_terms();
  term->set_coef(1.0);
  data::Determinant* det_hf = term->mutable_det();
  det_hf->mutable_up()->set_n_hf_elecs(n_up);
  det_hf->mutable_dn()->set_n_hf_elecs(n_dn);

  // Update HF and variational energy.
  energy_hf = energy_var = abstract_system->hamiltonian(det_hf, det_hf);

  // Clear connections.
  connections->clear();

  if (verbose) printf("HF det setup.\n");
};

void SolverImpl::variation(const double eps_var) {
  Timer* const timer = session->get_timer();

  // Contruct variational determinant hash set.
  data::Wavefunction* const wf = abstract_system->wf.get();
  const int n_start_dets = wf->terms_size();
  std::unordered_set<std::string> var_dets_set;
  std::vector<double> prev_coefs;
  var_dets_set.reserve(n_start_dets);
  prev_coefs.resize(n_start_dets, 0.0);
  for (const auto& term : wf->terms()) {
    var_dets_set.insert(term.det().SerializeAsString());
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
          var_dets_set.insert(std::move(det_code));
          data::Term* const term = wf->add_terms();
          term->set_coef(0.0);
          auto connected_det_copy = new data::Determinant(*connected_det);
          term->set_allocated_det(connected_det_copy);
          n_new_dets++;
        }
      };

  // Apply hamiltonian to a vec (with partially cached hamiltonian).
  const auto& apply_hamiltonian_func =
      std::bind(&SolverImpl::apply_hamiltonian, this, std::placeholders::_1);

  while (!converged) {
    timer->start("Variation: " + std::to_string(iteration));

    // Find new dets.
    n_new_dets = 0;
    const int n_old_dets = wf->terms_size();
    for (int i = 0; i < n_old_dets; i++) {
      const auto& term = wf->terms(i);
      const double coef = term.coef();
      if (std::abs(coef) <= std::abs(prev_coefs[i])) continue;
      abstract_system->find_connected_dets(
          &term.det(), eps_var / std::abs(term.coef()), connected_det_handler);
    }
    const int n_total_dets = wf->terms_size();
    if (verbose) {
      printf("New / total dets: %'d / %'d\n", n_new_dets, n_total_dets);
    }
    timer->checkpoint("found new dets");

    // Update prev coefs.
    prev_coefs.resize(n_total_dets);
    for (int i = 0; i < n_total_dets; i++) prev_coefs[i] = wf->terms(i).coef();

    // Diagonalize.
    connections->update();
    std::vector<double> diagonal(n_total_dets, 0.0);
    for (int i = 0; i < n_total_dets; i++) {
      const auto& det_i = wf->terms(i).det();
      diagonal[i] = abstract_system->hamiltonian(&det_i, &det_i);
    }
    const auto& diagonalization_result = DavidsonUtil::diagonalize(
        prev_coefs, diagonal, apply_hamiltonian_func, 10, verbose);
    const double energy_var_new = diagonalization_result.first;
    const auto& new_coefs = diagonalization_result.second;
    for (int i = 0; i < n_total_dets; i++) {
      wf->mutable_terms(i)->set_coef(new_coefs[i]);
    }

    // Determine convergence.
    if (std::abs(energy_var_new - energy_var) < 1.0e-6) {
      converged = true;
    }

    energy_var = energy_var_new;
    iteration++;

    if (verbose) print_var_result();

    timer->end();  // iteration.
  }
};

void SolverImpl::save_variation_result(const std::string& filename) {
  std::fstream var_file(
      filename, std::ios::out | std::ios::trunc | std::ios::binary);
  data::VariationResult res;
  res.set_energy_hf(energy_hf);
  res.set_energy_var(energy_var);
  res.set_allocated_wf(abstract_system->wf.release());
  res.SerializeToOstream(&var_file);
  var_file.close();
  abstract_system->wf.reset(res.release_wf());
  if (verbose) {
    printf("Saved to: %s\n", filename.c_str());
  }
  session->get_timer()->sleep(3);  // For data consistence on disk.
}

bool SolverImpl::load_variation_result(const std::string& filename) {
  std::fstream var_file(filename, std::ios::in | std::ios::binary);
  data::VariationResult res;
  if (!res.ParseFromIstream(&var_file)) {
    return false;
  }
  var_file.close();
  energy_hf = res.energy_hf();
  energy_var = res.energy_var();
  abstract_system->wf.reset(res.release_wf());
  if (verbose) {
    printf("Loaded from: %s\n", filename.c_str());
    print_var_result();
  }
  return true;
}

void SolverImpl::print_var_result() const {
  printf("Number of dets: %'d\n", abstract_system->wf->terms_size());
  printf("Variation energy: " ENERGY_FORMAT " Ha\n", energy_var);
  const double energy_corr = energy_var - energy_hf;
  printf("Correlation energy (var): " ENERGY_FORMAT " Ha\n", energy_corr);
}

std::vector<double> SolverImpl::apply_hamiltonian(
    const std::vector<double>& vec) {
  const int n_dets = vec.size();
  Parallel* const parallel = session->get_parallel();
  const int proc_id = parallel->get_proc_id();
  const int n_procs = parallel->get_n_procs();
  std::vector<double> res(n_dets, 0.0);

#pragma omp parallel for reduction(vec_double_plus : res) schedule(dynamic, 10)
  for (int i = proc_id; i < n_dets; i += n_procs) {
    const auto& conns = connections->get_connections(i);
    for (const auto conn : conns) {
      const int j = conn.first;
      const double H_ij = conn.second;
      res[i] += H_ij * vec[j];
      if (i != j) {
        res[j] += H_ij * vec[i];
      }
    }
  }

  parallel->reduce_to_sum(res);
  session->get_timer()->checkpoint("hamiltonian applied");

  return res;
};

void SolverImpl::perturbation(
    const int n_orbs_var,
    const double eps_var,
    const std::vector<int>& n_orbs_pts,
    const std::vector<double>& eps_pts) {
  // Clean variation variables.
  connections->clear();

  // Get Deterministic PT correction.
  const auto& energy_pts_dtm = get_energy_pts_dtm(n_orbs_pts);

  // Get Stochastic PT correction.
  const auto& energy_pts_stc =
      get_energy_pts_stc(n_orbs_pts, eps_pts, energy_pts_dtm);

  // Record results.
}

std::vector<double> SolverImpl::get_energy_pts_dtm(
    const std::vector<int>& n_orbs_pts) const {
  // Construct var dets set.
  std::unordered_set<std::string> var_dets_set;
  var_dets_set.reserve(abstract_system->wf->terms_size());
  for (const auto& term : abstract_system->wf->terms()) {
    var_dets_set.insert(term.det().SerializeAsString());
  }

  // Cache commonly used variables.
  const int n_n_orbs_pts = n_orbs_pts.size();
  const int n_var_dets = var_dets_set.size();
  Parallel* const parallel = session->get_parallel();
  Config* const config = session->get_config();
  Timer* const timer = session->get_timer();
  const size_t n_procs = parallel->get_n_procs();
  const int n_threads = parallel->get_n_threads();
  const size_t proc_id = parallel->get_proc_id();
  const size_t n_pt_batches_dtm = config->get_int("n_pt_batches_dtm");
  const double eps_pt_dtm = config->get_double("eps_pt_dtm");
  std::hash<std::string> string_hasher;

  // Construct partial sums store and pt results store.
  omp_hash_map<std::string, double> partial_sums;
  partial_sums.set_max_load_factor(1.618);
  std::vector<double> energy_pts_dtm(n_n_orbs_pts, 0.0);
  std::vector<unsigned long long> n_pt_dets_dtm(n_n_orbs_pts, 0);

  timer->start(str(boost::format("eps_pt_dtm: %#.4g") % eps_pt_dtm));
  // Process batch by batch to achieve a larger run in a constrained mem env.
  for (size_t b = 0; b < n_pt_batches_dtm; b++) {
    timer->start(
        str(boost::format("batch (%d/%d)") % (b + 1) % n_pt_batches_dtm));

    double target_progress = 0.25;
    timer->start("search");
#pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < n_var_dets; i++) {
      const auto& term = abstract_system->wf->terms(i);
      const auto& var_det = term.det();
      const double coef = term.coef();

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
        partial_sums.set(
            det_a_code, [&](double& value) { value += partial_sum_term; }, 0.0);
      };

      abstract_system->find_connected_dets(
          &var_det, eps_pt_dtm / std::abs(coef), pt_det_handler);

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
    timer->end();  // Search.

    timer->start("accumulate");

    // For reduction after each batch.
    std::vector<double> energy_pts_dtm_batch(n_n_orbs_pts, 0.0);

    // Avoid contructing new det every time.
    std::vector<data::Determinant> tmp_dets(n_threads);

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
      const double contribution = value * value / (energy_var - H_aa);
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
    for (int i = 0; i < n_n_orbs_pts; i++) {
      energy_pts_dtm[i] += energy_pts_dtm_batch[i];
    }

    timer->end();  // Accumulation.

    // Print batch result and estimated total correction.
    if (verbose) {
      printf("%20s", "CUM. energy PT:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf("%20.12f", energy_pts_dtm[i]);
      }
      printf("\n");
      printf("%20s", "EST. energy PT:");
      std::vector<double> energy_corrs_dtm(n_n_orbs_pts);
      for (int i = 0; i < n_n_orbs_pts; i++) {
        energy_corrs_dtm[i] = energy_pts_dtm[i] / (b + 1) * n_pt_batches_dtm;
        printf("%20.12f", energy_corrs_dtm[i]);
      }
      printf("\n");
      printf("%20s", "EST. COR. PT:");
      for (int i = 0; i < n_n_orbs_pts; i++) {
        printf("%20.12f", energy_corrs_dtm[i] + energy_var - energy_hf);
      }
      printf("\n");
    }

    partial_sums.clear();
    timer->end();  // Batch.
  }
  timer->end();

  parallel->reduce_to_sum(n_pt_dets_dtm);

  // Print total PT dets.
  if (verbose) {
    printf("%20s", "# DTM PT dets:");
    for (int i = 0; i < n_n_orbs_pts; i++) {
      printf("%'20llu", n_pt_dets_dtm[i]);
    }
    printf("\n");
  }

  return energy_pts_dtm;
}

std::vector<std::vector<double>> SolverImpl::get_energy_pts_stc(
    const std::vector<int>& n_orbs_pts,
    const std::vector<double>& eps_pts,
    const std::vector<double>& energy_pts_dtm) const {}

Solver* Injector::new_solver(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system) {
  return new SolverImpl(session, connections, abstract_system);
}
