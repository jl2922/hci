#ifndef SESSION_H_
#define SESSION_H_

#include "config.h"
#include "parallel.h"
#include "timer.h"

class Session {
 public:
  virtual ~Session() = default;

  virtual void init(int argc, char** argv) = 0;

  virtual Parallel* get_parallel() = 0;

  virtual Config* get_config() = 0;

  virtual Timer* get_timer() = 0;

  static Session* new_instance();
};

#endif