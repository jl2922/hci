#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

#include <vector>
#include "../abstract_system.h"
#include "../session.h"

class Connections {
 public:
  Connections(Session* const session, AbstractSystem* const abstract_system)
      : session(session), abstract_system(abstract_system) {}

  virtual ~Connections() = default;

  // Update by the new determinants.
  virtual void update() = 0;

  virtual void clear() = 0;

  virtual std::vector<std::pair<int, double>> get_connections(
      const int det_id) = 0;

 protected:
  Session* const session;

  AbstractSystem* const abstract_system;
};

#endif
