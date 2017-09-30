#include <clocale>
#include <memory>
#include "data.pb.h"
#include "injector.h"

int main(int argc, char** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  std::setlocale(LC_NUMERIC, "");

  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
