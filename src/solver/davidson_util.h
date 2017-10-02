#ifndef DAVIDSON_UTIL_H_
#define DAVIDSON_UTIL_H_

#include <functional>
#include <utility>
#include <vector>

class DavidsonUtil {
 public:
  static std::pair<double, std::vector<double>> diagonalize(
      const std::vector<double>& initial_vector,
      const std::vector<double>& diagonal,
      const std::function<std::vector<double>(const std::vector<double>&)>&
          apply_hamiltonian,
      const int max_iterations,
      const bool verbose);
};

#endif
