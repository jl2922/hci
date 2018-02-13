#include "chem_system.h"

#include <boost/functional/hash.hpp>
#include <cstdio>
#include <functional>
#include <numeric>
#include <vector>
#include "../injector.h"

#define ENERGY_FORMAT "%.12f"

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

  double energy_core;

  double energy_hf;

  std::unordered_map<size_t, double> f_integrals;

  std::unordered_map<size_t, double> g_integrals;

  std::vector<std::tuple<int, int, int, int, double>> raw_integrals;

  // std::unordered_map<std::array<int, 2>, double, boost::hash<std::array<int, 2>>> f_integral;

  // std::unordered_map<std::array<int, 4>, double, boost::hash<std::array<int, 4>>> g_integral;

  size_t combine_two(const size_t i, const size_t j);

  size_t get_integral_index(const size_t i, const size_t j, const size_t k, const size_t l);

  double integral_value(const int i, const int j, const int k, const int l);

  double integral_value(const int i, const int j);

  void load_fcidump();

  void reorder_orbitals();

  void obtain_hf();
};

ChemSystemImpl::ChemSystemImpl(Session* const session) : ChemSystem(session) {
  //...
}

void ChemSystemImpl::setup() {
  load_fcidump();

  reorder_orbitals();

  obtain_hf();
}

void ChemSystemImpl::load_fcidump() {
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
      energy_core = integral;
    } else if (r == s && s == 0) {
      f_integrals[combine_two(p - 1, q - 1)] = integral;
    } else {
      g_integrals[get_integral_index(p - 1, q - 1, r - 1, s - 1)] = integral;
    }
    raw_integrals.push_back(std::make_tuple(p, q, r, s, integral));
  }
}

void ChemSystemImpl::reorder_orbitals() {
  // Compute orbitals energy.
  std::vector<double> energy_orbs(n_orbs);
  for (int i = 0; i < n_orbs; i++) {
    energy_orbs[i] = integral_value(i, i);
    double energy_direct = 0;
    double energy_exchange = 0;
    for (int j = 0; j < n_up; j++) {
      if (hf_up[j] == i) continue;
      energy_exchange -= integral_value(i, hf_up[j], hf_up[j], i);
    }
    for (int j = 0; j < n_dn; j++) {
      if (hf_dn[j] == i) continue;
      energy_exchange -= integral_value(i, hf_dn[j], hf_dn[j], i);
    }
    for (int j = 0; j < n_up; j++) {
      const double integral = integral_value(i, i, hf_up[j], hf_up[j]);
      if (hf_up[j] == i) {
        energy_direct += integral;
      } else {
        energy_direct += 2 * integral;
      }
    }
    for (int j = 0; j < n_dn; j++) {
      const double integral = integral_value(i, i, hf_dn[j], hf_dn[j]);
      if (hf_dn[j] == i) {
        energy_direct += integral;
      } else {
        energy_direct += 2 * integral;
      }
    }
    energy_orbs[i] += 0.5 * (energy_direct + energy_exchange);
  }

  // Obtain new order.
  std::vector<int> orb_order(n_orbs);
  std::iota(orb_order.begin(), orb_order.end(), 0);
  std::sort(orb_order.begin(), orb_order.end(), [&](const int a, const int b) {
    return energy_orbs[a] < energy_orbs[b];
  });
  std::vector<int> orb_order_inv(n_orbs);
  for (int i = 0; i < n_orbs; i++) {
    orb_order_inv[orb_order[i]] = i;
  }
  printf("Orbitals reordered:\n");
  for (int i = 0; i < n_orbs; i++) {
    printf("new orb %d => old orb %d (E = %.12f)\n", i, orb_order[i], energy_orbs[orb_order[i]]);
  }

  // Update HF and integrals store.
  for (int i = 0; i < n_up; i++) {
    hf_up[i] = i;
  }
  for (int i = 0; i < n_dn; i++) {
    hf_dn[i] = i;
  }
  f_integrals.clear();
  g_integrals.clear();
  for (const auto& item : raw_integrals) {
    const int p = std::get<0>(item);
    const int q = std::get<1>(item);
    const int r = std::get<2>(item);
    const int s = std::get<3>(item);
    const double integral = std::get<4>(item);
    if (p == q && q == r && r == s && s == 0) {
      continue;
    } else if (r == s && s == 0) {
      f_integrals[combine_two(orb_order_inv[p - 1], orb_order_inv[q - 1])] = integral;
    } else {
      g_integrals[get_integral_index(
          orb_order_inv[p - 1], orb_order_inv[q - 1], orb_order_inv[r - 1], orb_order_inv[s - 1])] =
          integral;
    }
  }
  raw_integrals.clear();
}

void ChemSystemImpl::obtain_hf() {
  printf("Nuclear-Nuclear Energy: " ENERGY_FORMAT "\n", energy_core);
  double energy_one_body = 0;
  for (int i = 0; i < n_up; i++) {
    energy_one_body += f_integrals[combine_two(hf_up[i], hf_up[i])];
  }
  for (int i = 0; i < n_dn; i++) {
    energy_one_body += f_integrals[combine_two(hf_dn[i], hf_dn[i])];
  }
  printf("Energy One Body: " ENERGY_FORMAT "\n", energy_one_body);
  double energy_two_body = 0;
  double energy_direct = 0;
  double energy_exchange = 0;
  for (int i = 0; i < n_up; i++) {
    for (int j = i + 1; j < n_up; j++) {
      energy_direct += integral_value(hf_up[i], hf_up[i], hf_up[j], hf_up[j]);
      energy_exchange -= integral_value(hf_up[i], hf_up[j], hf_up[j], hf_up[i]);
    }
  }
  for (int i = 0; i < n_up; i++) {
    for (int j = 0; j < n_dn; j++) {
      energy_direct += integral_value(hf_up[i], hf_up[i], hf_dn[j], hf_dn[j]);
    }
  }
  for (int i = 0; i < n_dn; i++) {
    for (int j = i + 1; j < n_dn; j++) {
      energy_direct += integral_value(hf_dn[i], hf_dn[i], hf_dn[j], hf_dn[j]);
      energy_exchange -= integral_value(hf_dn[i], hf_dn[j], hf_dn[j], hf_dn[i]);
    }
  }
  printf("Energy Direct: " ENERGY_FORMAT "\n", energy_direct);
  printf("Energy Exchange: " ENERGY_FORMAT "\n", energy_exchange);
  energy_two_body = energy_direct + energy_exchange;
  printf("Energy Two Body: " ENERGY_FORMAT "\n", energy_two_body);
  energy_hf = energy_core + energy_one_body + energy_direct + energy_exchange;
  printf("Energy HF: " ENERGY_FORMAT "\n", energy_hf);
}

size_t ChemSystemImpl::combine_two(const size_t i, const size_t j) {
  if (i > j) {
    return (i * (i + 1)) / 2 + j;
  } else {
    return (j * (j + 1)) / 2 + i;
  }
}

size_t ChemSystemImpl::get_integral_index(
    const size_t i, const size_t j, const size_t k, const size_t l) {
  const size_t ij = combine_two(i, j);
  const size_t kl = combine_two(k, l);
  return combine_two(ij, kl);
}

double ChemSystemImpl::integral_value(const int i, const int j, const int k, const int l) {
  return g_integrals[get_integral_index(i, j, k, l)];
}

double ChemSystemImpl::integral_value(const int i, const int j) {
  return f_integrals[combine_two(i, j)];
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
