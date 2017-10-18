#ifndef ABSTRACT_SYSTEM_H_
#define ABSTRACT_SYSTEM_H_

#include <functional>
#include <memory>
#include "data.pb.h"
#include "session.h"

class AbstractSystem {
 public:
  AbstractSystem(Session* const session) : session(session) {}

  std::unique_ptr<data::Wavefunction> wf;

  virtual ~AbstractSystem() = default;

  virtual int get_n_orbitals() const = 0;

  virtual double hamiltonian(
      const data::Determinant* const det_pq,
      const data::Determinant* const det_rs) const = 0;

  virtual void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>&
          connected_det_handler) const = 0;

 protected:
  Session* const session;
};

#endif
