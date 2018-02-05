#pragma once

#include "../session.h"
#include "../solver/solver.h"
#include "chem_system.h"

class ChemController {
 public:
  ChemController(
      Session* const session,
      Solver* const solver,
      ChemSystem* const chem_system)
      : session(session), solver(solver), chem_system(chem_system) {}

  virtual ~ChemController() = default;

  virtual void run() = 0;

 protected:
  Session* const session;

  const std::unique_ptr<Solver> solver;

  const std::unique_ptr<ChemSystem> chem_system;
};
