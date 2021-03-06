#ifndef DAVIDSON_UTIL_H_
#define DAVIDSON_UTIL_H_

#include <functional>
#include <utility>
#include <vector>

class DavidsonUtil {
 public:
  // Returns lowest eigen value and eigen vector.
  static std::pair<double, std::vector<double>> diagonalize(
      const std::vector<double>& initial_vector,
      const std::vector<double>& diagonal,
      const std::function<std::vector<double>(
          const std::vector<double>&, const bool)>& apply_hamiltonian,
      const int max_iterations,
      const bool verbose);

 private:
  static void print_intermediate_result(
      const int iteration, const double lowest_eigenvalue);
};

#endif
