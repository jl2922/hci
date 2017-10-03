#include "connections.h"

#include "../injector.h"

class ConnectionsImpl : public Connections {
 public:
  ConnectionsImpl(
      Session* const session, AbstractSystem* const abstract_system);

  void update() override;

  std::vector<std::pair<int, double>> get_conections(const int i) override;

 private:
  Session* const session;

  AbstractSystem* const abstract_system;
};

ConnectionsImpl::ConnectionsImpl(
    Session* const session, AbstractSystem* const abstract_system)
    : session(session), abstract_system(abstract_system) {}

void ConnectionsImpl::update(){};

std::vector<std::pair<int, double>> ConnectionsImpl::get_conections(
    const int i) {
  std::vector<std::pair<int, double>> res;
  return res;
};

Connections* Injector::new_connections(
    Session* const session, AbstractSystem* const abstract_system) {
  return new ConnectionsImpl(session, abstract_system);
}
