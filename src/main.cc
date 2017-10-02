#include <clocale>
#include "data.pb.h"
#include "injector.h"

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  std::setlocale(LC_NUMERIC, "");

  std::unique_ptr<Session> session(Injector::new_session(argc, argv));

  Parallel* const parallel = session->get_parallel();
  if (parallel->is_master()) printf("\nHeat-bath Configuration Interaction\n");
  parallel->barrier();

  const auto& type = session->get_config()->get_string("type");
  if (type == "heg") {
    Injector::new_heg_controller(session.get())->run();
  } else {
    throw std::invalid_argument("System type '" + type + "' is not supported.");
  }

  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
