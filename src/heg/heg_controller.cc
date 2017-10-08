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

  std::vector<double> ruct_pts;

  std::vector<double> eps_pts;

  Session* const session;

  bool is_master = false;

  const std::unique_ptr<Solver> solver;

  const std::unique_ptr<HEGSystem> heg_system;

  void run_all_variations();

  void run_all_perturbations();
};

HEGControllerImpl::HEGControllerImpl(
    Session* const session, Solver* const solver, HEGSystem* const heg_system)
    : session(session), solver(solver), heg_system(heg_system) {
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
  for (int i = 0; i < n_rcut_vars; i++) {
    if (i > 0 && rcut_vars[i] == rcut_vars[i - 1]) continue;
    const double rcut_var = rcut_vars[i];
    timer->start(str(boost::format("rcut_var: %#.4g") % rcut_var));

    timer->start("setup");
    heg_system->setup(rcut_var);
    solver->setup_hf();
    timer->end();  // setup.

    for (int j = 0; j < n_eps_vars; j++) {
      if (j > 0) {
        if (eps_vars[j] > eps_vars[j - 1])
          throw std::invalid_argument("eps_var must be in decreasing order");
        if (eps_vars[j] == eps_vars[j - 1]) continue;
      }
      const double eps_var = eps_vars[j];
      timer->start(str(boost::format("eps_var: %#.4g") % eps_var));
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
  const double rcut_pt_ratio = config->get_double("rcut_pt_ratio");
  const double eps_pt_ratio = config->get_double("eps_pt_ratio");

  timer->start("perturbation");
  const int n_rcut_vars = rcut_vars.size();
  const int n_eps_vars = eps_vars.size();
  for (int i = n_rcut_vars - 1; i >= 0; i--) {
    const double rcut_var = rcut_vars[i];
    timer->start(str(boost::format("rcut_var: %#.4g") % rcut_var));
    const double rcut_pt = rcut_var * rcut_pt_ratio;
    if (is_master) printf("rcut_pt: %#.4g\n", rcut_pt);

    timer->start("setup");
    heg_system->setup(rcut_pt);
    timer->end();  // setup.

    for (int j = n_eps_vars - 1; j >= 0; j--) {
      const double eps_var = eps_vars[j];
      timer->start(str(boost::format("eps_var: %#.4g") % eps_var));
      const double eps_pt = eps_var * eps_pt_ratio;
      if (is_master) printf("eps_pt: %#.4g\n", eps_pt);

      const auto& filename =
          str(boost::format("var_%#.4g_%#.4g.dat") % eps_var % rcut_var);
      if (!solver->load_variation_result(filename)) {
        throw std::runtime_error("variational results missing.");
      }

      solver->perturbation(eps_pt);

      timer->end();  // eps_var.
    }

    timer->end();  // rcut_var.
  }
  timer->end();
}

HEGController* Injector::new_heg_controller(Session* const session) {
  HEGSystem* const heg_system = Injector::new_heg_system(session);
  Connections* const connections =
      Injector::new_connections(session, heg_system);
  Solver* const solver = Injector::new_solver(session, connections, heg_system);
  return Injector::new_heg_controller(session, solver, heg_system);
}

HEGController* Injector::new_heg_controller(
    Session* const session, Solver* const solver, HEGSystem* const heg_system) {
  return new HEGControllerImpl(session, solver, heg_system);
}
