#include <clocale>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "data.pb.h"
#include "injector.h"

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  std::setlocale(LC_NUMERIC, "");

  // Initialize parallel.
  Parallel* const parallel = Injector::new_parallel(argc, argv);

  // Initialize timer.
  Timer* const timer = Injector::new_timer(parallel);
  timer->init();

  // Initialize config.
  timer->start("Loading configuration");
  Config* const config = Injector::new_config("config.json", parallel);
  timer->end();

  // Initialize session.
  Session* const session = Injector::new_session(parallel, config, timer);

  const auto& type = config->get_string("type");
  if (type == "heg") {
    HEGSystem* const heg_system = Injector::new_heg_system(session);
    Connections* const connections =
        Injector::new_connections(session, heg_system);
    Solver* const solver =
        Injector::new_solver(session, connections, heg_system);
    HEGController* const heg_controller =
        Injector::new_heg_controller(session, solver, heg_system);
    heg_controller->run();
  } else if (type == "chem") {
    ChemSystem* const chem_system = Injector::new_chem_system(session);
    Connections* const connections =
        Injector::new_connections(session, chem_system);
    Solver* const solver =
        Injector::new_solver(session, connections, chem_system);
    ChemController* const chem_controller =
        Injector::new_chem_controller(session, solver, chem_system);
    chem_controller->run();
  } else {
    throw std::invalid_argument("System type '" + type + "' is not supported.");
  }

  delete session;

  google::protobuf::ShutdownProtobufLibrary();

  std::ofstream finish_log("finish.log", std::ios::out);

  return 0;
}
