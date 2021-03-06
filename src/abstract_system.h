#ifndef ABSTRACT_SYSTEM_H_
#define ABSTRACT_SYSTEM_H_

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "data.pb.h"
#include "session.h"

class AbstractSystem {
 public:
  AbstractSystem(Session* const session) : session(session) {}

  std::vector<std::string> dets;
  std::vector<double> coefs;

  virtual ~AbstractSystem() = default;

  virtual double hamiltonian(
      const data::Determinant* const det_pq,
      const data::Determinant* const det_rs) const = 0;

  virtual void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>&
          connected_det_handler) = 0;

 protected:
  Session* const session;
};

#endif
