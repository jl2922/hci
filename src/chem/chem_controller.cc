#include "chem_controller.h"

#include <boost/format.hpp>
#include <thread>
#include "../injector.h"

class ChemControllerImpl : public ChemController {
 public:
  ChemControllerImpl(
      Session* const session,
      Solver* const solver,
      ChemSystem* const chem_system);

  void run() override;

 private:
  std::vector<double> eps_vars;

  bool is_master = false;

  void run_all_variations();

  void run_all_perturbations();
};

ChemControllerImpl::ChemControllerImpl(
    Session* const session, Solver* const solver, ChemSystem* const chem_system)
    : ChemController(session, solver, chem_system) {
  is_master = session->get_parallel()->is_master();
}

void ChemControllerImpl::run() {
  run_all_variations();

  if (session->get_config()->get_bool("variation_only")) return;

  std::this_thread::sleep_for(std::chrono::seconds(1));

  run_all_perturbations();
}

void ChemControllerImpl::run_all_variations() {
  Timer* const timer = session->get_timer();
  Config* const config = session->get_config();
  eps_vars = config->get_double_array("eps_vars");
  const int n_eps_vars = eps_vars.size();

  timer->start("variation");

  // Check eps in decreasing order.
  for (int i = 1; i < n_eps_vars; i++) {
    assert(eps_vars[i] < eps_vars[i - 1]);
  }

  timer->start("setup");
  chem_system->setup();
  solver->setup_hf();
  timer->end();  // setup.

  for (int j = 0; j < n_eps_vars; j++) {
    const double eps_var = eps_vars[j];
    timer->start(str(boost::format("eps_var %#.4g") % eps_var));
    // const auto& filename =
    //     str(boost::format("var_%#.4g_%#.4g.dat") % eps_var % rcut_var);
    // if (!solver->load_variation_result(filename)) {
    //   solver->variation(eps_var);
    //   if (is_master) solver->save_variation_result(filename);
    // }
    timer->end();  // eps_var.
  }

  timer->end();  // variation.
}

void ChemControllerImpl::run_all_perturbations() {}

ChemController* Injector::new_chem_controller(
    Session* const session,
    Solver* const solver,
    ChemSystem* const chem_system) {
  return new ChemControllerImpl(session, solver, chem_system);
}
