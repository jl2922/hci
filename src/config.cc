#include "config.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <string>
#include "injector.h"

class ConfigImpl : public Config {
 public:
  ConfigImpl(const std::string& filename, Parallel* const parallel);

  std::string get_string(const std::string& key) const override {
    return config_tree.get<std::string>(key);
  }

  bool get_bool(const std::string& key) const override {
    return config_tree.get<bool>(key);
  }

  int get_int(const std::string& key) const override {
    return config_tree.get<int>(key);
  }

  double get_double(const std::string& key) const override {
    return config_tree.get<double>(key);
  }

  std::vector<double> get_double_array(const std::string& key) const override;

 private:
  boost::property_tree::ptree config_tree;
};

ConfigImpl::ConfigImpl(const std::string& filename, Parallel* const parallel) {
  boost::property_tree::read_json(filename, config_tree);
  if (parallel->is_master()) {
    printf("Configuration:\n");
    boost::property_tree::json_parser::write_json(std::cout, config_tree);
  }
  parallel->barrier();
}

std::vector<double> ConfigImpl::get_double_array(const std::string& key) const {
  std::vector<double> res;
  for (auto& item : config_tree.get_child(key)) {
    res.push_back(item.second.get_value<double>());
  }
  return res;
}

Config* Injector::new_config(
    const std::string& filename, Parallel* const parallel) {
  return new ConfigImpl(filename, parallel);
}