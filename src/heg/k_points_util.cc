#include "k_points_util.h"

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <exception>
#include <unordered_set>
#include "../array_math.h"

int KPointsUtil::get_n_k_points(const double rcut) {
  int count = 0;
  const int n_max = std::floor(static_cast<float>(rcut));
  double rcut_square = rcut * rcut;
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > rcut_square) continue;
        count++;
      }
    }
  }
  return count;
}

std::vector<std::array<int8_t, 3>> KPointsUtil::generate_k_points(
    const double rcut) {
  std::vector<std::array<int8_t, 3>> k_points;
  const int n_max = std::floor(static_cast<float>(rcut));
  assert(n_max < 127);
  double rcut_square = rcut * rcut;
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > rcut_square) continue;
        k_points.push_back(std::array<int8_t, 3>({static_cast<int8_t>(i),
                                                  static_cast<int8_t>(j),
                                                  static_cast<int8_t>(k)}));
      }
    }
  }
  std::stable_sort(
      k_points.begin(),
      k_points.end(),
      [](const auto& a, const auto& b) -> bool {
        return squared_norm(a) < squared_norm(b);
      });
  return k_points;
}

std::vector<std::array<int8_t, 3>> KPointsUtil::get_k_diffs(
    const std::vector<std::array<int8_t, 3>>& k_points) {
  // Generate all possible differences between two different k points.
  std::unordered_set<std::array<int8_t, 3>, boost::hash<std::array<int8_t, 3>>>
      k_diffs_set;
  std::vector<std::array<int8_t, 3>> k_diffs;
  const int n_k_points = k_points.size();
  for (int p = 0; p < n_k_points; p++) {
    for (int q = 0; q < n_k_points; q++) {
      if (p == q) continue;
      const auto& diff_pq = k_points[q] - k_points[p];
      if (k_diffs_set.count(diff_pq) == 1) continue;
      k_diffs.push_back(diff_pq);
      k_diffs_set.insert(diff_pq);
    }
  }

  // Sort k_diffs into ascending order so that later sorting hci queue will be
  // faster.
  std::stable_sort(
      k_diffs.begin(), k_diffs.end(), [](const auto& a, const auto& b) -> bool {
        return squared_norm(a) < squared_norm(b);
      });

  return k_diffs;
}