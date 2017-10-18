#ifndef HEG_SYSTEM_H_
#define HEG_SYSTEM_H_

#include "../abstract_system.h"

class HEGSystem : public AbstractSystem {
 public:
  HEGSystem(Session* const session) : AbstractSystem(session) {}

  virtual ~HEGSystem() = default;

  // Generate k points, hci array, and HF for the specified rcut.
  virtual void setup(const double rcut) = 0;

  virtual int get_n_orbitals(const double rcut) const = 0;
};

#endif
