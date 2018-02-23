#include "chem_system.h"

#include <boost/functional/hash.hpp>
#include <cstdio>
#include <functional>
#include <numeric>
#include <vector>
#include "../injector.h"
#include "product_table.h"

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

  std::string point_group;

  int n_group_elements;

  int ms2;

  std::vector<int> hf_up;

  std::vector<int> hf_dn;

  std::vector<int> orb_syms;

  double energy_core;

  double energy_hf;

  bool is_dih;

  std::unordered_map<size_t, double> f_integrals;

  std::unordered_map<size_t, double> g_integrals;

  std::vector<std::vector<int>> product_table;

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

  void init_point_group();

  void init_product_table(const int** product_table_array);

  void convert_to_adams_indices();
};

ChemSystemImpl::ChemSystemImpl(Session* const session) : ChemSystem(session) {
  //...
}

void ChemSystemImpl::setup() {
  load_fcidump();

  reorder_orbitals();

  obtain_hf();

  init_point_group();
}

void ChemSystemImpl::load_fcidump() {
  FILE* fcidump = fopen("FCIDUMP", "r");
  if (!fcidump) throw new std::runtime_error("FCIDUMP not found");
  char buf[80];
  fscanf(fcidump, "%*s %*s %d", &n_orbs);
  printf("N_ORBS: %d\n", n_orbs);
  orb_syms.resize(n_orbs);
  int n_elecs;
  fscanf(fcidump, "%*s %d", &n_elecs);
  printf("N_ELECS: %d\n", n_elecs);
  n_up = n_dn = n_elecs / 2;
  fscanf(fcidump, "%*s %d %*s", &ms2);
  printf("MS2: %d\n", ms2);
  fscanf(fcidump, "  ORBSYM=");
  for (int i = 0; i < n_orbs; i++) {
    int sym;
    fscanf(fcidump, "%d,", &sym);
    orb_syms[i] = sym;
  }
  fscanf(fcidump, "%*s");
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
  std::vector<int> orb_syms_new(n_orbs);
  for (int i = 0; i < n_orbs; i++) {
    orb_syms_new[i] = orb_syms[orb_order[i]];
  }
  orb_syms = std::move(orb_syms_new);
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

void ChemSystemImpl::init_point_group() {
  is_dih = false;
  Config* const config = session->get_config();
  point_group = config->get_string("point_group");
  printf("PG: %s\n", point_group.c_str());
  if (point_group == "c1") {
    n_group_elements = 1;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::C1[0], ProductTable::C1[0] + n_group_elements);
    }
  } else if (point_group == "cs") {
    n_group_elements = 2;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::CS[0], ProductTable::CS[0] + n_group_elements);
    }
  } else if (point_group == "ci") {
    n_group_elements = 2;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::CI[0], ProductTable::CI[0] + n_group_elements);
    }
  } else if (point_group == "c2v") {
    n_group_elements = 4;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::C2V[0], ProductTable::C2V[0] + n_group_elements);
    }
  } else if (point_group == "c2h") {
    n_group_elements = 4;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::C2H[0], ProductTable::C2H[0] + n_group_elements);
    }
  } else if (point_group == "d2h") {
    n_group_elements = 8;
    product_table.resize(n_group_elements);
    for (int i = 0; i < n_group_elements; i++) {
      product_table[i].assign(ProductTable::D2H[0], ProductTable::D2H[0] + n_group_elements);
    }
  } else if (point_group == "dih") {
    n_group_elements = 18;
    convert_to_adams_indices();
  } else {
    throw std::invalid_argument("The point group is not implemented");
  }
}

void ChemSystemImpl::convert_to_adams_indices() {
  bool is_adams_indices = true;
  for (int i = 0; i < n_orbs; i++) {
    if (orb_syms[i] < 0) {
      is_adams_indices = false;
      break;
    }
  }
  if (!is_adams_indices) {
    std::unordered_map<int, int> sandeep_to_adams_map;
    for (int i = 0; i < n_orbs; i++) {
      sandeep_to_adams_map[ProductTable::SandeepDIHIndices[i]] = i + 1;
    }
    for (int i = 0; i < n_orbs; i++) {
      printf("%d => %d\n ", orb_syms[i],  sandeep_to_adams_map[orb_syms[i]]);
      orb_syms[i] = sandeep_to_adams_map[orb_syms[i]];
    }
  }
}

// case ('c1')
//    n_group_elements = 1
//    allocate(product_table_elems(n_group_elements, n_group_elements))
//    product_table_elems(1,:) = (/1 /)
// case ('cs')
//    n_group_elements = 2
//    allocate(product_table_elems(n_group_elements, n_group_elements))
//    product_table_elems(1,:) = (/1, 2 /)
//    product_table_elems(2,:) = (/2, 1 /)
// case ('ci')

//    ! Sandeep's indices:  1, 2, 5, 6,-5,-6, 7, 8,-7, -8,  9, 10, -9,-10, 11, 12,-11,-12
//    ! Adam's indices:     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
//    ! angular momentum:   0, 0, 1, 1,-1,-1, 2, 2,-2, -2,  3,  3, -3, -3,  4,  4, -4, -4
//    ! g/u (g=0,u=1):      0, 1, 0, 1, 0, 1, 0, 1, 0,  1,  0,  1,  0,  1,  0,  1,  0,  1

//    if (minval(orbital_symmetries)<0) then ! there are negative indices, so assume that the
//    indices given are Sandeep's. This converts them to my indices.
//      write(6,*)
//      write(6,'(''Replacing orbital sym 1, 2, 5, 6,-5,-6, 7, 8,-7, -8,  9, 10, -9,-10, 11,
//      12,-11,-12  by'')') write(6,'(''                      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
//      13, 14, 15, 16, 17, 18'')') write(6,'(''angular momentum:     0, 0, 1, 1,-1,-1, 2, 2,-2, -2,
//      3,  3, -3, -3,  4,  4, -4, -4'')') write(6,'(''g/u (g=0,u=1):        0, 1, 0, 1, 0, 1, 0, 1,
//      0,  1,  0,  1,  0,  1,  0,  1,  0,  1'',/)') do i=1,norb
//        if (orbital_symmetries(i).ne.1.and.orbital_symmetries(i).ne.2) then
//          old = orbital_symmetries(i)
//          a = abs(old)/2
//          b = (abs(old)+1)/2
//          new = a+3*b-8
//          if (old<0)  new=new+2
//          orbital_symmetries(i) = new
//        endif
//      enddo
//    endif

//    call get_lz(maxval(orbital_symmetries),lz,gu)
//    !Max number of group elements is based off of angular momentum selection rules
//    !We must take 3 * maximum angular momentum, as that is the maximum angular momentum
//    !possible, and 4 * that, as there are 4 indices per angular momentum,
//    !and an additional 2 for the 0 lz state.

//    n_group_elements=12*abs(lz)+2

//    d_infinity_h = 1

//    ! Following product_table is implemented in the code in a way that allows for arbitrarily
//    ! large L_z:

//   !n_group_elements = 18
//   !allocate(product_table_elems(n_group_elements, n_group_elements))
//   !product_table_elems(1,:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18/)
//   !product_table_elems(2,:) = (/ 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
//   !product_table_elems(3,:) = (/ 3, 4, 7, 8, 1, 2,11,12, 5, 6,15,16, 9,10,19,20,13,14/)
//   !product_table_elems(4,:) = (/ 4, 3, 8, 7, 2, 1,12,11, 6, 5,16,15,10, 9,20,19,14,13/)
//   !product_table_elems(5,:) = (/ 5, 6, 1, 2, 9,10, 3, 4,13,14, 7, 8,17,18,11,12,21,22/)
//   !product_table_elems(6,:) = (/ 6, 5, 2, 1,10, 9, 4, 3,14,13, 8, 7,18,17,12,11,22,21/)
//   !product_table_elems(7,:) = (/ 7, 8,11,12, 3, 4,15,16, 1, 2,19,20, 5, 6,23,24, 9,10/)
//   !product_table_elems(8,:) = (/ 8, 7,12,11, 4, 3,16,15, 2, 1,20,19, 6, 5,24,23,10, 9/)
//   !product_table_elems(9,:) = (/ 9,10, 5, 6,13,14, 1, 2,17,18, 3, 4,21,22, 7, 8,25,26/)
//   !product_table_elems(10,:)= (/10, 9, 6, 5,14,13, 2, 1,18,17, 4, 3,22,21, 8, 7,26,25/)
//   !product_table_elems(11,:)= (/11,12,15,16, 7, 8,19,20, 3, 4,23,24, 1, 2,27,28, 5, 6/)
//   !product_table_elems(12,:)= (/12,11,16,15, 8, 7,20,19, 4, 3,24,23, 2, 1,28,27, 6, 5/)
//   !product_table_elems(13,:)= (/13,14, 9,10,17,18, 5, 6,21,22, 1, 2,25,26, 3, 4,29,30/)
//   !product_table_elems(14,:)= (/14,13,10, 9,18,17, 6, 5,22,21, 2, 1,26,25, 4, 3,30,29/)
//   !product_table_elems(15,:)= (/15,16,19,20,11,12,23,24, 7, 8,27,28, 3, 4,31,32, 1, 2/)
//   !product_table_elems(16,:)= (/16,15,20,19,12,11,24,23, 8, 7,28,27, 4, 3,32,31, 2, 1/)
//   !product_table_elems(17,:)= (/17,18,13,14,21,22, 9,10,25,26, 5, 6,29,30, 1, 2,33,34/)
//   !product_table_elems(18,:)= (/18,17,14,13,22,21,10, 9,26,25, 6, 5,30,29, 2, 1,34,33/)

// case default
//    write(6,*) "The point group ", point_group, " is not implemented"
//    stop "The point group is not implemented"
// end select

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
