#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

#include <vector>
#include "../data.pb.h"

class Connections {
 public:
  virtual ~Connections() = default;

  // Update by the new determinants.
  virtual void update() = 0;

  virtual void clear() = 0;

  virtual std::vector<std::pair<int, double>> get_connections(const int i) = 0;

  virtual std::vector<std::pair<int, double>> get_connections(
      const data::Determinant& det, const int i, const double eps) = 0;
};

#endif
