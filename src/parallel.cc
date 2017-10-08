#include "parallel.h"

#include <boost/mpi.hpp>
#include "injector.h"
#include "omp.h"
#include "vector_plus.h"

class ParallelImpl : public Parallel {
 public:
  ParallelImpl(int argc, char** argv);

  bool is_master() const override { return proc_id == 0; }

  int get_proc_id() const override { return proc_id; }

  int get_n_procs() const override { return n_procs; }

  int get_n_threads() const override { return n_threads; }

  void barrier() override {
    fflush(stdout);
    world.barrier();
  }

  void reduce_to_sum(unsigned long long& value) {
    unsigned long long local_value = value;
    boost::mpi::all_reduce(
        world, local_value, value, std::plus<unsigned long long>());
  }

  void reduce_to_sum(double& value) {
    double local_value = value;
    boost::mpi::all_reduce(world, local_value, value, std::plus<double>());
  }

  void reduce_to_sum(std::vector<unsigned long long>& value) {
    std::vector<unsigned long long> local_value = value;
    boost::mpi::all_reduce(
        world, local_value, value, VectorPlus<unsigned long long>());
  }

  void reduce_to_sum(std::vector<double>& value) {
    std::vector<double> local_value = value;
    boost::mpi::all_reduce(world, local_value, value, VectorPlus<double>());
  }

 private:
  int proc_id = 0;

  int n_procs = 0;

  int n_threads = 0;

  std::unique_ptr<boost::mpi::environment> env;

  boost::mpi::communicator world;
};

ParallelImpl::ParallelImpl(int argc, char** argv) {
  env.reset(new boost::mpi::environment(argc, argv));
  proc_id = world.rank();
  n_procs = world.size();
  n_threads = omp_get_max_threads();
  omp_set_nested(1);

  if (is_master()) printf("Infrastructure:\n");
  barrier();
  printf(
      "Proc #%d (%d threads): %s\n",
      proc_id,
      n_threads,
      env->processor_name().c_str());
  barrier();
}

Parallel* Injector::new_parallel(int argc, char** argv) {
  return new ParallelImpl(argc, argv);
}
