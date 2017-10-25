#include "heg_controller.h"

#include <boost/format.hpp>
#include <exception>
#include <memory>
#include "../injector.h"

class HEGControllerImpl : public HEGController {
 public:
  HEGControllerImpl(
      Session* const session,
      Solver* const solver,
      HEGSystem* const heg_system);

  void run() override;

 private:
  std::vector<double> rcut_vars;

  std::vector<double> eps_vars;

  bool is_master = false;

  void run_all_variations();

  void run_all_perturbations();
};

HEGControllerImpl::HEGControllerImpl(
    Session* const session, Solver* const solver, HEGSystem* const heg_system)
    : HEGController(session, solver, heg_system) {
  is_master = session->get_parallel()->is_master();
}

void HEGControllerImpl::run() {
  run_all_variations();

  if (session->get_config()->get_bool("variation_only")) return;

  run_all_perturbations();
}

void HEGControllerImpl::run_all_variations() {
  Timer* const timer = session->get_timer();
  Config* const config = session->get_config();
  rcut_vars = config->get_double_array("rcut_vars");
  eps_vars = config->get_double_array("eps_vars");

  timer->start("variation");
  const int n_rcut_vars = rcut_vars.size();
  const int n_eps_vars = eps_vars.size();

  // Check eps in decreasing order.
  for (int i = 1; i < n_eps_vars; i++) {
    if (eps_vars[i] >= eps_vars[i - 1])
      throw std::invalid_argument("eps_var must be in decreasing order");
  }

  for (int i = 0; i < n_rcut_vars; i++) {
    const double rcut_var = rcut_vars[i];
    timer->start(str(boost::format("rcut_var %#.4g") % rcut_var));

    timer->start("setup");
    heg_system->setup(rcut_var);
    solver->setup_hf();
    timer->end();  // setup.

    for (int j = 0; j < n_eps_vars; j++) {
      const double eps_var = eps_vars[j];
      timer->start(str(boost::format("eps_var %#.4g") % eps_var));
      const auto& filename =
          str(boost::format("var_%#.4g_%#.4g.dat") % eps_var % rcut_var);
      if (!solver->load_variation_result(filename)) {
        solver->variation(eps_var);
        if (is_master) solver->save_variation_result(filename);
      }
      timer->end();  // eps_var.
    }

    timer->end();  // rcut_var.
  }

  timer->end();  // variation.
}

void HEGControllerImpl::run_all_perturbations() {
  Timer* const timer = session->get_timer();
  Config* const config = session->get_config();
  const auto& rcut_pts = config->get_double_array("rcut_pts");
  const auto& eps_pts = config->get_double_array("eps_pts");
  const int n_rcut_vars = rcut_vars.size();
  const int n_eps_vars = eps_vars.size();
  const int n_rcut_pts = rcut_pts.size();
  const int n_eps_pts = eps_pts.size();
  std::vector<int> n_orbs_pts;

  for (int i = 0; i < n_rcut_pts; i++) {
    const double rcut_pt = rcut_pts[i];
    n_orbs_pts.push_back(heg_system->get_n_orbitals(rcut_pt));
  }

  timer->start("perturbation");

  for (int i = 1; i < n_eps_pts; i++) {
    if (eps_pts[i] >= eps_pts[i - 1])
      throw std::invalid_argument("eps_pts must be in decreasing order.");
  }

  const double max_rcut_pt =
      *std::max_element(rcut_pts.begin(), rcut_pts.end());

  if (is_master) {
    printf("Maximum PT rcut: %#.4g\n", max_rcut_pt);
    printf("Minimum PT epsilon: %#.4g\n", eps_pts[n_eps_pts - 1]);
  }

  timer->start("setup");
  heg_system->setup(max_rcut_pt);
  timer->end();

  for (int i = n_rcut_vars - 1; i >= 0; i--) {
    const double rcut_var = rcut_vars[i];
    const int n_orbs_var = heg_system->get_n_orbitals(rcut_var);
    timer->start(str(boost::format("rcut_var %#.4g") % rcut_var));

    for (int j = n_eps_vars - 1; j >= 0; j--) {
      const double eps_var = eps_vars[j];
      timer->start(str(boost::format("eps_var %#.4g") % eps_var));
      const auto& filename =
          str(boost::format("var_%#.4g_%#.4g.dat") % eps_var % rcut_var);
      if (!solver->load_variation_result(filename)) {
        throw std::runtime_error("variational results missing.");
      }
      solver->perturbation(n_orbs_var, eps_var, n_orbs_pts, eps_pts);
      timer->end();  // eps_var.
    }

    timer->end();  // rcut_var.
  }

  timer->end();
}

HEGController* Injector::new_heg_controller(
    Session* const session, Solver* const solver, HEGSystem* const heg_system) {
  return new HEGControllerImpl(session, solver, heg_system);
}
