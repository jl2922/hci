#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <algorithm>
#include <functional>
#include <vector>

class Parallel {
 public:
  Parallel(int argc, char** argv) {
    (void)argc;
    (void)argv;
  };

  virtual ~Parallel() = default;

  virtual bool is_master() const = 0;

  virtual int get_proc_id() const = 0;

  virtual int get_n_procs() const = 0;

  virtual int get_n_threads() const = 0;

  virtual void barrier() = 0;

  virtual void reduce_to_sum(unsigned long long& value) = 0;

  virtual void reduce_to_sum(double& value) = 0;

  virtual void reduce_to_sum(std::vector<unsigned long long>& value) = 0;

  virtual void reduce_to_sum(std::vector<double>& value) = 0;

  virtual void reduce_to_sum(std::vector<long double>& value) = 0;
};

#endif
