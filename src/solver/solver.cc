#include "solver.h"

#include <cmath>
#include <string>
#include <unordered_set>
#include <utility>
#include "../injector.h"
#include "connections.h"
#include "davidson_util.h"

#define ENERGY_FORMAT "%.12f"

class SolverImpl : public Solver {
 public:
  SolverImpl(
      Session* const session,
      Connections* const connections,
      AbstractSystem* const abstrct_system);

  void setup_hf() override;

  void variation(const double eps_var) override;

  std::vector<double> apply_hamiltonian(
      const std::vector<double>& vec) override;

  void perturbation(const double eps_pt) override;

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
}

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

    if (verbose) {
      const double energy_corr = energy_var - energy_hf;
      printf("Variation energy: " ENERGY_FORMAT " Ha\n", energy_var_new);
      printf("Correlation energy: " ENERGY_FORMAT " Ha\n", energy_corr);
    }

    timer->end();  // iteration.
  }
};

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

void SolverImpl::perturbation(const double eps_pt){};

Solver* Injector::new_solver(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system) {
  return new SolverImpl(session, connections, abstract_system);
}
