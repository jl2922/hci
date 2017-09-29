#ifndef ABSTRACT_SYSTEM_H_
#define ABSTRACT_SYSTEM_H_

#include "data.pb.h"

class AbstractSystem {
 public:
  virtual ~AbstractSystem() = default;

  virtual data::Wavefunction* get_wf() = 0;

  virtual double hamiltonian(
      const data::Determinant* const det_pq, const data::Determinant* const det_rs) const = 0;

  virtual void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>& connected_det_handler) const = 0;
};

#endif