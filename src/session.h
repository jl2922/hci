#ifndef SESSION_H_
#define SESSION_H_

#include "config.h"
#include "parallel.h"
#include "timer.h"

class Session {
 public:
  Session(Parallel* const parallel, Config* const config, Timer* const timer)
      : parallel(parallel), config(config), timer(timer) {}

  virtual ~Session() = default;

  virtual Parallel* get_parallel() = 0;

  virtual Config* get_config() = 0;

  virtual Timer* get_timer() = 0;

 protected:
  const std::unique_ptr<Parallel> parallel;

  const std::unique_ptr<Config> config;

  const std::unique_ptr<Timer> timer;
};

#endif
