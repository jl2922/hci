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

  void perturbation(const double eps_pt) override;

  void perturbation(
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts) override;

 private:
  int n_up = 0;

  int n_dn = 0;

  double energy_hf = 0.0;

  double energy_var = 0.0;

  double energy_pt = 0.0;

  bool verbose = false;

  Session* const session;

  const std::unique_ptr<Connections> connections;

  AbstractSystem* const abstract_system;

  std::ofstream pt_result;

  void print_var_result() const;
};

SolverImpl::SolverImpl(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system)
    : session(session),
      connections(connections),
      abstract_system(abstract_system) {
  verbose = session->get_parallel()->is_master();
  Config* const config = session->get_config();
  n_up = config->get_int("n_up");
  n_dn = config->get_int("n_dn");
  if (!config->get_bool("variation_only")) {
    pt_result.open("pt_result.csv", std::ios::out | std::ios::trunc);
    pt_result << "n_orbs_pt,eps_pt,n_pt_dets,energy_corr" << std::endl;
  }
}

SolverImpl::~SolverImpl() { pt_result.close(); }

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
  data::SpinDeterminant* up = det_hf->mutable_up();
  up->set_n_hf_elecs(n_up);
  data::SpinDeterminant* dn = det_hf->mutable_dn();
  dn->set_n_hf_elecs(n_dn);

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
    std::string det_code;
    term.det().SerializeToString(&det_code);
    var_dets_set.insert(std::move(det_code));
  }

  // Variation iterations.
  bool converged = false;
  int iteration = 0;
  int n_new_dets = 0;
  const auto& connected_det_handler =
      [&](const data::Determinant* const connected_det) {
        std::string det_code;
        connected_det->SerializeToString(&det_code);
        if (var_dets_set.count(det_code) == 0) {
          var_dets_set.insert(std::move(det_code));
          data::Term* const term = wf->add_terms();
          term->set_coef(0.0);
          auto connected_det_copy = new data::Determinant(*connected_det);
          term->set_allocated_det(connected_det_copy);
          n_new_dets++;
        }
      };
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
    print_var_result();
    printf("Saved to: %s\n", filename.c_str());
  }
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

void SolverImpl::perturbation(const double eps_pt) {
  std::vector<int> n_orbs_pts;
  std::vector<double> eps_pts;
  n_orbs_pts.push_back(abstract_system->get_n_orbitals());
  eps_pts.push_back(eps_pt);
  perturbation(n_orbs_pts, eps_pts);
};

void SolverImpl::perturbation(
    const std::vector<int>& n_orbs_pts, const std::vector<double>& eps_pts) {
  // Clean variation variables.
  connections->clear();

  // Construct var dets hash set.
  std::unordered_set<std::string> var_dets_set;
  var_dets_set.reserve(abstract_system->wf->terms_size());
  for (const auto& term : abstract_system->wf->terms()) {
    var_dets_set.insert(term.det().SerializeAsString());
  }

  // Cache variables.
  const int n_n_orbs_pts = n_orbs_pts.size();
  const int n_eps_pts = eps_pts.size();
  const double min_eps_pt = eps_pts.back();
  const int n_var_dets = var_dets_set.size();
  Parallel* const parallel = session->get_parallel();
  Timer* const timer = session->get_timer();
  const int n_procs = parallel->get_n_procs();
  const int n_threads = parallel->get_n_threads();
  const int proc_id = parallel->get_proc_id();
  std::hash<std::string> string_hasher;

  // Setup partial sum hash map.
  std::vector<omp_hash_map<std::string, double>> partial_sums(n_eps_pts);
  for (int i = 0; i < n_eps_pts; i++) {
    partial_sums[i].set_max_load_factor(1.618);
  }

  // Search pt dets.
  double target_progress = 0.25;
  timer->start("search");
#pragma omp parallel for schedule(static, 1)
  for (int i = 0; i < n_var_dets; i++) {
    const auto& term = abstract_system->wf->terms(i);
    const auto& var_det = term.det();
    const double coef = term.coef();

    const auto& pt_det_handler = [&](const auto& det_a) {
      const auto& det_a_code = det_a->SerializeAsString();
      const int det_a_hash = string_hasher(det_a_code);
      if (var_dets_set.count(det_a_code) == 1) return;
      if (det_a_hash % n_procs != proc_id) return;
      const double H_ai = abstract_system->hamiltonian(&var_det, det_a);
      const double partial_sum_term = H_ai * coef;
      for (int j = 0; j < n_eps_pts; j++) {
        if (std::abs(partial_sum_term) >= eps_pts[j]) {
          partial_sums[j].set(
              det_a_code,
              [&](double& value) { value += partial_sum_term; },
              0.0);
          break;
        }
      }
    };

    abstract_system->find_connected_dets(
        &var_det, min_eps_pt / std::abs(coef), pt_det_handler);

    // Report progress.
    const double current_progress = i * 100.0 / n_var_dets;
    if (target_progress <= current_progress) {
      const int thread_id = omp_get_thread_num();
      if (thread_id == 0) {
        std::string event =
            str(boost::format("Progress: %.2f %%") % target_progress);
        timer->checkpoint(event);
        if (verbose) {
          unsigned long long master_n_dets_total = 0;
          for (int j = 0; j < n_eps_pts; j++) {
            master_n_dets_total += partial_sums[j].get_n_keys();
          }
          printf("Master node PT dets: %'llu\n", master_n_dets_total);
        }
        target_progress *= 2.0;
      }
    }
  }
  timer->end();

  // Accumulate partial results.
  timer->start("accumulation");
  std::vector<std::vector<unsigned long long>> n_pt_dets(n_eps_pts);
  std::vector<std::vector<double>> energy_pts(n_eps_pts);
  for (int i = 0; i < n_eps_pts; i++) {
    n_pt_dets[i].resize(n_n_orbs_pts, 0);
    energy_pts[i].resize(n_n_orbs_pts, 0.0);
  }
  std::vector<data::Determinant> tmp_dets(n_threads);
  for (int i = 0; i < n_eps_pts; i++) {
    partial_sums[i].apply([&](const std::string& det_code, const double value) {
      const int thread_id = omp_get_thread_num();
      auto& det = tmp_dets[thread_id];
      det.ParseFromString(det_code);
      const int n_orbs_used = 0;
      const double H_aa = abstract_system->hamiltonian(&det, &det);
      double partial_sum_value = value;
      for (int j = i; j < n_eps_pts; j++) {
        if (j > i) {
          const double partial_sum_j =
              partial_sums[j].get_copy_or_default(det_code, 0.0);
          partial_sums[j].unset(det_code);
          partial_sum_value += partial_sum_j;
        }
        const double contribution =
            partial_sum_value * partial_sum_value / (energy_var - H_aa);
        for (int k = 0; k < n_n_orbs_pts; k++) {
          if (n_orbs_used > n_orbs_pts[k]) continue;
#pragma omp atomic
          n_pt_dets[j][k]++;
#pragma omp atomic
          energy_pts[j][k] += contribution;
        }
      }
    });
  }

  // Report results.
  for (int i = 0; i < n_eps_pts; i++) {
    parallel->reduce_to_sum(n_pt_dets[i]);
    parallel->reduce_to_sum(energy_pts[i]);
    const double cur_eps_pt = eps_pts[i];

    if (verbose) {
      for (int j = 0; j < n_n_orbs_pts; j++) {
        if (i != 0 || j != 0) printf("\n");
        const int cur_n_orbs_pt = n_orbs_pts[j];
        const double cur_energy_pt = energy_pts[i][j];
        const unsigned long long cur_n_pt_dets = n_pt_dets[i][j];

        printf("PT eps: %#.4g\n", cur_eps_pt);
        printf("PT n orbs: %'d\n", cur_n_orbs_pt);
        printf("Perturbation energy: %#.15g Ha\n", cur_energy_pt);
        const double energy_corr = cur_energy_pt + energy_var - energy_hf;
        printf("Correlation energy (pt): %#.15g Ha\n", energy_corr);
        printf("Number of perturbation dets: %'llu\n", cur_n_pt_dets);
        // pt_result << str(boost::format("%d, %#.4g, %llu, %.12f") % n_pt_orbs
        // %
        //                  eps_pt % n_pt_dets % energy_corr)
        //           << std::endl;
      }
    }
  }

  timer->end();
}

Solver* Injector::new_solver(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system) {
  return new SolverImpl(session, connections, abstract_system);
}
