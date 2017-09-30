#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <cstddef>
#include <vector>

class Parallel {
 public:
  virtual ~Parallel() = default;

  virtual bool is_master() const = 0;

  virtual int get_proc_id() const = 0;

  virtual int get_n_procs() const = 0;

  virtual int get_n_threads() const = 0;

  virtual void barrier() = 0;

  virtual void reduce_to_sum(unsigned long long& value) = 0;

  virtual void reduce_to_sum(double& value) = 0;

  virtual void reduce_to_sum(std::vector<double>& value) = 0;

  static Parallel* new_instance(int argc, char** argv);
};

#endif
