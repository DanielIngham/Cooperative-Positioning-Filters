/**
 * @file measurements.hpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */

#include "CL/sensors/meas_set.hpp"

#include "UtiasMrclam/agents/Robot.hpp"

#include <vector>

namespace CL::sensors {

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

  MeasSet *const find(double time);

private:
  std::vector<MeasSet> measurements_;
};
} // namespace CL::sensors
