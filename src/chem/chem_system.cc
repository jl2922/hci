#include "chem_system.h"

#include <boost/functional/hash.hpp>
#include <cstdio>
#include "../injector.h"

class ChemSystemImpl : public ChemSystem {
 public:
  ChemSystemImpl(Session* const session);

  void setup() override;

  double hamiltonian(
      const data::Determinant* const det_pq, const data::Determinant* const det_rs) const override;

  void find_connected_dets(
      const data::Determinant* const det,
      const double eps,
      const std::function<void(const data::Determinant* const)>& connected_det_handler) override;

 private:
  int n_orbs;

  int n_up;

  int n_dn;

  int ms2;

  std::vector<int> hf_up;

  std::vector<int> hf_dn;

  double E_nn;

  std::unordered_map<std::array<int, 2>, double, boost::hash<std::array<int, 2>>> f_int;

  std::unordered_map<std::array<int, 4>, double, boost::hash<std::array<int, 4>>> g_int;
};

ChemSystemImpl::ChemSystemImpl(Session* const session) : ChemSystem(session) {
  //...
}

void ChemSystemImpl::setup() {
  FILE* fcidump = fopen("FCIDUMP", "r");
  if (!fcidump) throw new std::runtime_error("FCIDUMP not found");
  char buf[80];
  fscanf(fcidump, "%*s %*s %d", &n_orbs);
  printf("N_ORBS: %d\n", n_orbs);
  int n_elecs;
  fscanf(fcidump, "%*s %d", &n_elecs);
  printf("N_ELECS: %d\n", n_elecs);
  n_up = n_dn = n_elecs / 2;
  fscanf(fcidump, "%*s %d %*s", &ms2);
  printf("MS2: %d\n", ms2);
  fscanf(fcidump, "%*s %*s");
  hf_up.resize(n_up);
  hf_dn.resize(n_dn);
  printf("HF UP: ");
  for (int i = 0; i < n_up; i++) {
    fscanf(fcidump, "%d", &hf_up[i]);
    hf_up[i]--;
    printf("%d ", hf_up[i]);
  }
  fscanf(fcidump, "%*s");
  printf("\nHF DN: ");
  for (int i = 0; i < n_dn; i++) {
    fscanf(fcidump, "%d", &hf_dn[i]);
    hf_dn[i]--;
    printf("%d ", hf_dn[i]);
  }
  fscanf(fcidump, "%*s %*s");
  printf("\n");
  int p, q, r, s;
  double integral;
  while (fscanf(fcidump, "%lf %d %d %d %d", &integral, &p, &q, &r, &s) != EOF) {
    if (p == q && q == r && r == s && s == 0) {
      E_nn = integral;
    } else if (r == s && s == 0) {
      f_int[std::array<int, 2>{p - 1, q - 1}] = integral;
    } else {
      g_int[std::array<int, 4>{p - 1, q - 1, r - 1, s - 1}] = integral;
    }
  }
}

double ChemSystemImpl::hamiltonian(
    const data::Determinant* const det_pq, const data::Determinant* const det_rs) const {
  return 0.0;
}

void ChemSystemImpl::find_connected_dets(
    const data::Determinant* const det,
    const double eps,
    const std::function<void(const data::Determinant* const)>& connected_det_handler) {}

ChemSystem* Injector::new_chem_system(Session* const session) {
  return new ChemSystemImpl(session);
}
