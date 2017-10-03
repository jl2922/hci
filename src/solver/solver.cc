#include "solver.h"

#include "../injector.h"

class SolverImpl : public Solver {
 public:
  SolverImpl(
      Session* const session,
      Connections* const connections,
      AbstractSystem* const abstrct_system);

  void setup_hf() override;

  void variation(const double eps_var) override;

  void perturbation(const double eps_pt) override;

 private:
  Session* const session;

  const std::unique_ptr<Connections> connections;

  AbstractSystem* const abstract_system;
};

SolverImpl::SolverImpl(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstrct_system)
    : session(session),
      connections(connections),
      abstract_system(abstract_system) {}

void SolverImpl::setup_hf(){};

void SolverImpl::variation(const double eps_var){};

void SolverImpl::perturbation(const double eps_pt){};

Solver* Injector::new_solver(
    Session* const session,
    Connections* const connections,
    AbstractSystem* const abstract_system) {
  return new SolverImpl(session, connections, abstract_system);
}
