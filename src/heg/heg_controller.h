#ifndef HEG_CONTROLLER_H_
#define HEG_CONTROLLER_H_

#include "../session.h"
#include "../solver/solver.h"
#include "heg_system.h"

class HEGController {
 public:
  HEGController(
      Session* const session, Solver* const solver, HEGSystem* const heg_system)
      : session(session), solver(solver), heg_system(heg_system) {}

  virtual ~HEGController() = default;

  virtual void run() = 0;

 protected:
  Session* const session;

  const std::unique_ptr<Solver> solver;

  const std::unique_ptr<HEGSystem> heg_system;
};

#endif
