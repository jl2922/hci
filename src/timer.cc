#include "timer.h"

#include <chrono>
#include <cstdio>
#include <ctime>
#include <thread>
#include "injector.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

class TimerImpl : public Timer {
 public:
  TimerImpl(Parallel* const parallel) : Timer(parallel) {
    verbose = parallel->is_master();
  }

  void init() override;

  void start(const std::string& event) override;

  void checkpoint(const std::string& msg) override;

  void end() override;

 private:
  bool verbose = false;

  std::chrono::high_resolution_clock::time_point init_time;

  std::chrono::high_resolution_clock::time_point prev_time;

  std::vector<
      std::pair<std::string, std::chrono::high_resolution_clock::time_point>>
      start_times;

  void barrier() { parallel->barrier(); }

  void print_event_path() const;

  void print_time() const;

  double get_duration(
      const std::chrono::high_resolution_clock::time_point start,
      const std::chrono::high_resolution_clock::time_point end) const {
    return (std::chrono::duration_cast<std::chrono::duration<double>>(
                end - start))
        .count();
  }
};

void TimerImpl::init() {
  barrier();
  const auto& now = std::chrono::high_resolution_clock::now();
  init_time = prev_time = now;
  if (verbose) {
    const time_t start_time = std::chrono::system_clock::to_time_t(now);
    printf("\nStart time: %s", asctime(localtime(&start_time)));
    printf("Format: " ANSI_COLOR_YELLOW "[DIFF/SECTION/TOTAL]" ANSI_COLOR_RESET
           "\n");
  }
  barrier();
}

void TimerImpl::start(const std::string& event) {
  barrier();
  const auto& now = std::chrono::high_resolution_clock::now();
  start_times.push_back(std::make_pair(event, now));
  if (verbose) {
    printf("\n" ANSI_COLOR_GREEN "[START] " ANSI_COLOR_RESET);
    print_event_path();
    printf(" ");
    print_time();
    printf("\n");
  }
  prev_time = now;
  barrier();
}

void TimerImpl::checkpoint(const std::string& msg) {
  barrier();
  const auto& now = std::chrono::high_resolution_clock::now();
  if (verbose) {
    printf(ANSI_COLOR_GREEN "[-CHK-] " ANSI_COLOR_RESET "%s ", msg.c_str());
    print_time();
    printf("\n");
  }
  prev_time = now;
  barrier();
}

void TimerImpl::end() {
  barrier();
  const auto& now = std::chrono::high_resolution_clock::now();
  if (verbose) {
    printf(ANSI_COLOR_GREEN "[=END=] " ANSI_COLOR_RESET);
    print_event_path();
    printf(" ");
    print_time();
    printf("\n");
  }
  start_times.pop_back();
  prev_time = now;
  barrier();
}

void TimerImpl::print_event_path() const {
  for (size_t i = 0; i < start_times.size() - 1; i++) {
    printf("%s >> ", start_times[i].first.c_str());
  }
  printf("%s", start_times.back().first.c_str());
}

void TimerImpl::print_time() const {
  const auto& now = std::chrono::high_resolution_clock::now();
  const auto& event_start_time = start_times.back().second;
  printf(
      ANSI_COLOR_YELLOW "[%.3f/%.3f/%.3f]" ANSI_COLOR_RESET,
      get_duration(prev_time, now),
      get_duration(event_start_time, now),
      get_duration(init_time, now));
}

Timer* Injector::new_timer(Parallel* const parallel) {
  return new TimerImpl(parallel);
}
