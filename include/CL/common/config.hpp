/**
 * @file config.hpp
 */
#pragma once

#include <cstddef>
#include <string>

namespace CL {

/**
 * Datastructure housing all inference related configuration parameters.
 */
class Config {
public:
  Config() = default;
  Config(Config &&) = default;
  Config(const Config &) = default;
  Config &operator=(Config &&) = default;
  Config &operator=(const Config &) = default;
  ~Config() = default;

  /**
   * Extracts a given configuration yaml file.
   * @param config_file Path to the configuration file.
   */
  Config(const std::string &config_file);

  /**
   * Houses the parameters that the inference class would use to configure the
   * behaviour of the robot agents in the VANET.
   */
  struct Robots {
    /** The number of cooperative robots. */
    size_t cooperative{1};
    /** The number of adversarial robots. */
    size_t adversarial{};
    /** The number of faulty cooperative robots. */
    size_t faulty{};
  } robots; /**< Configuration for the robots in the VANET. */

  /**
   * Houses the parameters that the inference class would use to configure the
   * behaviour of the landmark agents in the VANET.
   */
  struct Landmarks {
    /** The number of cooperative robots. */
    size_t cooperative{15};
    /** The number of faulty cooperative robots. */
    size_t adversarial{};
  } landmarks; /**< Configuration for the landmarks in the VANET. */

private:
};
} // namespace CL
