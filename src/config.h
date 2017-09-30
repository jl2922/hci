#ifndef CONFIG_H_
#define CONFIG_H_

#include <memory>
#include <string>
#include <vector>
#include "parallel.h"

class Config {
 public:
  virtual ~Config() = default;

  virtual std::string get_string(const std::string& key) const = 0;

  virtual bool get_bool(const std::string& key) const = 0;

  virtual int get_int(const std::string& key) const = 0;

  virtual double get_double(const std::string& key) const = 0;

  virtual std::vector<double> get_double_array(const std::string& key) const = 0;
};

#endif
