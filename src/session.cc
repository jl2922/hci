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
};

SessionImpl::SessionImpl(
    Parallel* const parallel, Config* const config, Timer* const timer)
    : Session(parallel, config, timer) {}

Session* Injector::new_session(
    Parallel* const parallel, Config* const config, Timer* const timer) {
  return new SessionImpl(parallel, config, timer);
}
