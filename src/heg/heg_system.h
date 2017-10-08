#ifndef HEG_SYSTEM_H_
#define HEG_SYSTEM_H_

#include "../abstract_system.h"

class HEGSystem : public AbstractSystem {
 public:
  virtual ~HEGSystem() = default;

  virtual void setup(const double rcut) = 0;

  virtual int get_n_orbitals(const double rcut) const = 0;
};

#endif
