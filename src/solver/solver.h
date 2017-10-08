#ifndef SOLVER_H_
#define SOLVER_H_

#include "../abstract_system.h"
#include "../session.h"

class Solver {
 public:
  virtual ~Solver() = default;

  virtual void setup_hf() = 0;

  virtual void variation(const double eps_var) = 0;

  virtual void save_variation_result(const std::string& filename) = 0;

  virtual bool load_variation_result(const std::string& filename) = 0;

  virtual std::vector<double> apply_hamiltonian(
      const std::vector<double>& vec) = 0;

  virtual void perturbation(
      const int n_orbs_var,  // For output.
      const double eps_var,  // For output.
      const std::vector<int>& n_orbs_pts,
      const std::vector<double>& eps_pts) = 0;
};

#endif
