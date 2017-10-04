#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

#include <vector>

class Connections {
 public:
  virtual ~Connections() = default;

  virtual void update() = 0;

  virtual void clear() = 0;

  virtual std::vector<std::pair<int, double>> get_connections(const int i) = 0;
};

#endif
