#ifndef TIMER_H_
#define TIMER_H_

#include "parallel.h"

#include <string>

class Timer {
 public:
  virtual ~Timer() = default;

  virtual void init() = 0;

  virtual void start(const std::string& event) = 0;

  virtual void checkpoint(const std::string& msg) = 0;

  virtual void end() = 0;

  static Timer* new_instance(Parallel* const parallel);
};

#endif