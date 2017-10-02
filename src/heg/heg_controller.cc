#include "heg_controller.h"

#include "../injector.h"

class HEGControllerImpl : public HEGController {
 public:
  HEGControllerImpl(
      Session* const session,
      Solver* const solver,
      HEGSystem* const heg_system);

  void run() override;

 private:
  Session* const session;

  Solver* const solver;

  HEGSystem* const heg_system;

  void variation();

  void perturbation();
};

HEGControllerImpl::HEGControllerImpl(
    Session* const session, Solver* const solver, HEGSystem* const heg_system)
    : session(session), solver(solver), heg_system(heg_system) {}

void HEGControllerImpl::run() {
  variation();

  perturbation();
}

void HEGControllerImpl::variation() {
  Timer* const timer = session->get_timer();

  timer->start("variation");
  timer->end();
}

void HEGControllerImpl::perturbation() {
  Timer* const timer = session->get_timer();

  timer->start("perturbation");
  timer->end();
}

HEGController* Injector::new_heg_controller(Session* const session) {
  return Injector::new_heg_controller(session, nullptr, nullptr);
}

HEGController* Injector::new_heg_controller(
    Session* const session, Solver* const solver, HEGSystem* const heg_system) {
  return new HEGControllerImpl(session, solver, heg_system);
}
