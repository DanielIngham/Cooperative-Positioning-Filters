#include "CL/common/config.hpp"

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/node/parse.h>
#include <yaml-cpp/yaml.h>

namespace CL {

Config::Config(const std::string &config_file) {
  if (!std::filesystem::exists(config_file)) {
    throw std::invalid_argument(
        "Unable to find specified configuration file: " + config_file);
  }

  YAML::Node config{YAML::LoadFile(config_file)};

  robots = {
      .faulty = config["robots"]["faulty"].as<size_t>(),
  };

  landmarks = {
      .adversarial = config["landmarks"]["adversarial"].as<size_t>(),
  };
}
} // namespace CL
