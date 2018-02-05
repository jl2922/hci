#include "chem_system.h"

#include "../injector.h"

class ChemSystemImpl : public ChemSystem {
 public:
  ChemSystemImpl(Session* const session);

  void setup() override;

  double hamiltonian(
      const data::Determinant* const det_pq,
      const data::Determinant* const det_rs) const override;

  void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>&
          connected_det_handler) override;
};

ChemSystemImpl::ChemSystemImpl(Session* const session) : ChemSystem(session) {
  //...
}

void ChemSystemImpl::setup() {}

double ChemSystemImpl::hamiltonian(
    const data::Determinant* const det_pq,
    const data::Determinant* const det_rs) const {
  return 0.0;
}

void ChemSystemImpl::find_connected_dets(
    const data::Determinant* const det,
    const double eps,
    const std::function<void(const data::Determinant* const)>&
        connected_det_handler) {}

ChemSystem* Injector::new_chem_system(Session* const session) {
  return new ChemSystemImpl(session);
}
