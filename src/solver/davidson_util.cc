#include "davidson_util.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

std::pair<double, std::vector<double>> DavidsonUtil::diagonalize(
    const std::vector<double>& initial_vector,
    const std::vector<double>& diagonal,
    const std::function<std::vector<double>(const std::vector<double>&)>&
        apply_hamiltonian,
    const int max_iterations_in,
    const bool verbose) {
  const double TOLERANCE = 1.0e-10;
  const double EPSILON = 1.0e-10;
  const int n_dets = initial_vector.size();
  double lowest_eigenvalue = 0.0;
  std::vector<double> lowest_eigenvector;

  if (n_dets == 1) {
    lowest_eigenvalue = diagonal[0];
    lowest_eigenvector.resize(1);
    lowest_eigenvector[0] = 1.0;
    return std::make_pair(lowest_eigenvalue, std::move(lowest_eigenvector));
  }

  const int max_iterations = std::min(n_dets, max_iterations_in);
  double lowest_eigenvalue_prev = 0.0;
  double residual_norm = 0.0;

  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(n_dets, max_iterations);
  for (int i = 0; i < n_dets; i++) v(i, 0) = initial_vector[i];
  v.col(0).normalize();

  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(n_dets, max_iterations);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(n_dets);
  Eigen::VectorXd Hw = Eigen::VectorXd::Zero(n_dets);
  Eigen::MatrixXd h_krylov =
      Eigen::MatrixXd::Zero(max_iterations, max_iterations);
  Eigen::MatrixXd h_overwrite;
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(max_iterations);
  int len_work = 3 * max_iterations - 1;
  Eigen::VectorXd work(len_work);
  bool converged = false;
  std::vector<double> tmp_v(n_dets);

  // Get diagonal elements.
  Eigen::VectorXd diag_elems(n_dets);
  for (int i = 0; i < n_dets; i++) diag_elems[i] = diagonal[i];

  // First iteration.
  for (int i = 0; i < n_dets; i++) tmp_v[i] = v(i, 0);
  const auto& tmp_Hv = apply_hamiltonian(tmp_v);
  for (int i = 0; i < n_dets; i++) Hv(i, 0) = tmp_Hv[i];
  lowest_eigenvalue = v.col(0).dot(Hv.col(0));
  h_krylov(0, 0) = lowest_eigenvalue;
  w = v.col(0);
  Hw = Hv.col(0);
  if (verbose) print_intermediate_result(1, lowest_eigenvalue);

  for (int it = 1; it < max_iterations; it++) {
    // Compute residual.
    for (int j = 0; j < n_dets; j++) {
      const double diff_to_diag = lowest_eigenvalue - diag_elems[j];
      if (std::abs(diff_to_diag) < EPSILON) {
        v(j, it) = -1.0;
      } else {
        v(j, it) = (Hw(j, 0) - lowest_eigenvalue * w(j, 0)) / diff_to_diag;
      }
    }

    // If residual is small, converge.
    residual_norm = v.col(it).norm();
    if (residual_norm < TOLERANCE) converged = true;

    // Orthogonalize and normalize.
    for (int i = 0; i < it; i++) {
      double norm = v.col(it).dot(v.col(i));
      v.col(it) -= norm * v.col(i);
    }
    v.col(it).normalize();

    // Apply H once.
    for (int i = 0; i < n_dets; i++) tmp_v[i] = v(i, it);
    const auto& tmp_Hv2 = apply_hamiltonian(tmp_v);
    for (int i = 0; i < n_dets; i++) Hv(i, it) = tmp_Hv2[i];

    // Construct subspace matrix.
    for (int i = 0; i <= it; i++) {
      h_krylov(i, it) = v.col(i).dot(Hv.col(it));
      h_krylov(it, i) = h_krylov(i, it);
    }

    // Diagonalize subspace matrix.
    len_work = 3 * it + 2;
    h_overwrite = h_krylov.leftCols(it + 1).topRows(it + 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(
        h_krylov.leftCols(it + 1).topRows(it + 1));
    const auto& eigenvalues = eigenSolver.eigenvalues();
    const auto& eigenvectors = eigenSolver.eigenvectors();
    lowest_eigenvalue = eigenvalues[0];
    int lowest_id = 0;
    for (int i = 1; i <= it; i++) {
      if (eigenvalues[i] < lowest_eigenvalue) {
        lowest_eigenvalue = eigenvalues[i];
        lowest_id = i;
      }
    }
    w = v.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
    Hw = Hv.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);

    if (std::abs(lowest_eigenvalue - lowest_eigenvalue_prev) < TOLERANCE) {
      converged = true;
      break;
    } else {
      lowest_eigenvalue_prev = lowest_eigenvalue;
      if (verbose) print_intermediate_result(it, lowest_eigenvalue);
    }

    if (converged) break;
  }

  lowest_eigenvector.resize(n_dets);
  for (int i = 0; i < n_dets; i++) lowest_eigenvector[i] = w(i);

  return std::make_pair(lowest_eigenvalue, std::move(lowest_eigenvector));
}

void DavidsonUtil::print_intermediate_result(
    const int iteration, const double lowest_eigenvalue) {
  printf(
      "Davidson Iteration #%d. Eigenvalue: %#.12f\n",
      iteration,
      lowest_eigenvalue);
}
