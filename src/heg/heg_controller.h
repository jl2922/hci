#ifndef HEG_CONTROLLER_H_
#define HEG_CONTROLLER_H_

#include "../session.h"
#include "../solver/solver.h"
#include "heg_system.h"

class HEGController {
 public:
  virtual ~HEGController() = default;

  // Read config.
  // For each rcut_var in rcut_vars:
  //   Setup heg system.
  //   Run Variation for each eps_var in eps_vars.
  // For each rcut_var and eps_var:
  //   Run perturbation.
  virtual void run() = 0;

  static HEGController* new_instance(Session* const session);
};

#endif
