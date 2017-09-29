#ifndef HEG_CONTROLLER_H_
#define HEG_CONTROLLER_H_

#include "../session.h"

class HEGController {
 public:
  virtual ~HEGController() = default;

  virtual void run() = 0;

  static HEGController* new_instance(Session* const session);
};

#endif