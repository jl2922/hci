#ifndef ARRAY_MATH_H_
#define ARRAY_MATH_H_

#include <array>

template <class T, size_t N>
std::array<T, N> operator+(
    const std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  std::array<T, N> res;
  for (size_t i = 0; i < N; i++) res[i] = lhs[i] + rhs[i];
  return res;
}

template <class T, size_t N>
std::array<T, N> operator-(
    const std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  std::array<T, N> res;
  for (size_t i = 0; i < N; i++) res[i] = lhs[i] - rhs[i];
  return res;
}

template <class T, size_t N>
std::array<double, N> operator*(const std::array<T, N> &lhs, const double rhs) {
  std::array<double, N> res;
  for (size_t i = 0; i < N; i++) res[i] = lhs[i] * rhs;
  return res;
}

template <class T, size_t N>
double squared_norm(const std::array<T, N> &arr) {
  double res = 1.0 * arr[0] * arr[0];  // Prevent overflow.
  for (size_t i = 1; i < N; i++) res += 1.0 * arr[i] * arr[i];
  return res;
}

template <class T, size_t N>
double norm(const std::array<T, N> &arr) {
  return sqrt(squared_norm(arr));
}

template <class T1, class T2, size_t N>
bool operator==(const std::array<T1, N> &lhs, const T2 &rhs) {
  for (size_t i = 0; i < N; i++) {
    if (lhs[i] != rhs) return false;
  }
  return true;
}

template <class T1, class T2, size_t N>
bool operator!=(const std::array<T1, N> &lhs, const T2 &rhs) {
  return !(lhs == rhs);
}

#endif
