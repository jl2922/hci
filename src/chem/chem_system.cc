#include "chem_system.h"

#include <cstdio>
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

 private:
  int n_orbs;

  int n_up;

  int n_dn;
};

ChemSystemImpl::ChemSystemImpl(Session* const session) : ChemSystem(session) {
  //...
}

void ChemSystemImpl::setup() {
  FILE* fcidump = fopen("FCIDUMP", "r");
  if (!fcidump) throw new std::runtime_error("FCIDUMP not found");
  char buf[80];
  fscanf(fcidump, "%*s %*s %d", &n_orbs);
  printf("%d\n", n_orbs);
}

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
