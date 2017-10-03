#include "davidson_util.h"

std::pair<double, std::vector<double>> DavidsonUtil::diagonalize(
    const std::vector<double>& initial_vector,
    const std::vector<double>& diagonal,
    const std::function<std::vector<double>(const std::vector<double>&)>&
        apply_hamiltonian,
    const int max_iterations,
    const bool verbose) {
  const int n_dets = initial_vector.size();
  std::pair<double, std::vector<double>> res;
  res.first = diagonal[0];
  res.second.resize(n_dets);
  return res;
}
