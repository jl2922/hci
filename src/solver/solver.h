#ifndef SOLVER_H_
#define SOLVER_H_

#include "../abstract_system.h"
#include "../session.h"
#include "connections.h"

class Solver {
 public:
  Solver(
      Session* const session,
      Connections* const connections,
      AbstractSystem* const abstract_system)
      : session(session),
        connections(connections),
        abstract_system(abstract_system) {
    parallel = session->get_parallel();
    config = session->get_config();
    timer = session->get_timer();
    n_procs = parallel->get_n_procs();
    n_threads = parallel->get_n_threads();
    proc_id = parallel->get_proc_id();
  }

  virtual ~Solver() = default;

  virtual void setup_hf() = 0;

  virtual void variation(const double eps_var) = 0;

  virtual void save_variation_result(const std::string& filename) = 0;

  virtual bool load_variation_result(const std::string& filename) = 0;

  virtual std::vector<double> apply_hamiltonian(
      const std::vector<double>& vec, const bool first_iteration) = 0;

  virtual void perturbation(
      const int n_orbs_var,
      const double eps_var,
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts) = 0;

 protected:
  Session* const session;

  Parallel* parallel;

  Config* config;

  Timer* timer;

  size_t n_procs;

  int n_threads;

  size_t proc_id;

  const std::unique_ptr<Connections> connections;

  AbstractSystem* const abstract_system;
};

#endif
