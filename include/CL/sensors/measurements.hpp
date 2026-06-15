/**
 * @file measurements.hpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */

#include "UtiasMrclam/agents/Robot.hpp"
#include <vector>

class Measurements {
public:
  Measurements() = delete;
  Measurements(Measurements &&) = default;
  Measurements(const Measurements &) = default;
  Measurements &operator=(Measurements &&) = default;
  Measurements &operator=(const Measurements &) = default;
  ~Measurements() = default;

  Measurements(
      std::vector<utias::mrclam::Robot::Measurement> const &measurements);

private:
};
