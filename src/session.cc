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

Session* Injector::new_session(int argc, char** argv) {
  // Initialize parallel.
  Parallel* const parallel = Injector::new_parallel(argc, argv);

  // Initialize timer.
  Timer* const timer = Injector::new_timer(parallel);
  timer->init();

  // Initialize config.
  timer->start("Loading configuration");
  Config* const config = Injector::new_config("config.json", parallel);
  timer->end();

  return Injector::new_session(parallel, config, timer);
}

Session* Injector::new_session(
    Parallel* const parallel, Config* const config, Timer* const timer) {
  return new SessionImpl(parallel, config, timer);
}
