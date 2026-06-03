#include "CL/common/config.hpp"

#include <cstddef>
#include <filesystem>
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
      .cooperative = config["robots"]["cooperative"].as<size_t>(),
      .faulty = config["robots"]["cooperative"].as<size_t>(),
  };

  landmarks = {
      .cooperative = config["landmarks"]["cooperative"].as<size_t>(),
      .adversarial = config["landmarks"]["adversarial"].as<size_t>(),
  };
}
} // namespace CL
