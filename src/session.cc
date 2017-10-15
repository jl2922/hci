#include "session.h"

#include <memory>
#include "injector.h"

class SessionImpl : public Session {
 public:
  SessionImpl(
      Parallel* const parallel, Config* const config, Timer* const timer);

  Parallel* get_parallel() override { return parallel.get(); };

  Config* get_config() override { return config.get(); };

  Timer* get_timer() override { return timer.get(); };

 private:
  const std::unique_ptr<Parallel> parallel;

  const std::unique_ptr<Config> config;

  const std::unique_ptr<Timer> timer;
};

SessionImpl::SessionImpl(
    Parallel* const parallel, Config* const config, Timer* const timer)
    : parallel(parallel), config(config), timer(timer) {}

Session* Injector::new_session(
    Parallel* const parallel, Config* const config, Timer* const timer) {
  return new SessionImpl(parallel, config, timer);
}
