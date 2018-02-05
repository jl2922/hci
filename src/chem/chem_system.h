#pragma once

#include "../abstract_system.h"

class ChemSystem : public AbstractSystem {
 public:
  ChemSystem(Session* const session) : AbstractSystem(session) {}

  virtual ~ChemSystem() = default;

  // Read FCIDUMP, generate hci arrays, and calculate HF energy.
  virtual void setup() = 0;
};
