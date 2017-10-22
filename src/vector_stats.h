#ifndef VECTOR_STATS_H_
#define VECTOR_STATS_H_

#include <vector>

template <class T>
double get_avg(const std::vector<T>& vec) {
  double sum_elems = 0.0;
  const int n_elems = vec.size();
  for (int i = 0; i < n_elems; i++) sum_elems += vec[i];
  return sum_elems / n_elems;
}

template <class T>
double get_stdev(const std::vector<T>& vec) {
  double sum_vars = 0.0;
  const double avg = get_avg(vec);
  const int n_elems = vec.size();
  for (int i = 0; i < n_elems; i++) {
    const double diff = avg - vec[i];
    sum_vars += diff * diff;
  }
  return sqrt(sum_vars / (n_elems - 1));
}

#endif
