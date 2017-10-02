#ifndef INJECTOR_H_
#define INJECTOR_H_

#include "abstract_system.h"
#include "config.h"
#include "heg/heg_controller.h"
#include "heg/heg_system.h"
#include "parallel.h"
#include "session.h"
#include "solver/connections.h"
#include "solver/solver.h"
#include "timer.h"

// Dependency injection.
class Injector {
 public:
  static Session* new_session(int argc, char** argv);

  static Session* new_session(
      Parallel* const parallel, Config* const config, Timer* const timer);

  static Parallel* new_parallel(int argc, char** argv);

  static Config* new_config(
      const std::string& filename, Parallel* const parallel);

  static Timer* new_timer(Parallel* const parallel);

  static Solver* new_solver(
      Session* const session,
      Connections* const connections,
      AbstractSystem* const system);

  static Connections* new_connections(
      Session* const session, AbstractSystem* const system);

  static HEGController* new_heg_controller(Session* const session);

  static HEGController* new_heg_controller(
      Session* const session,
      Solver* const solver,
      HEGSystem* const heg_system);

  static HEGSystem* new_heg_system(Session* const session);
};

#endif
