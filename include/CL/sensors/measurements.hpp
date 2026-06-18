/**
 * @file measurements.hpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */

#include "CL/common/types.hpp"
#include "CL/sensors/meas_set.hpp"

#include "UtiasMrclam/agents/Robot.hpp"

#include <vector>

namespace CL::sensors {

/**
 * Houses all the measurement data for a given agent.
 */
class Measurements {
public:
  Measurements() = delete;
  Measurements(Measurements &&) = default;
  Measurements(const Measurements &) = default;
  Measurements &operator=(Measurements &&) = default;
  Measurements &operator=(const Measurements &) = default;
  ~Measurements() = default;

  /**
   * Construct a list of measurements observed by the agent.
   * @param measurements List of measurements for the robot provided by the data
   * handler.
   * @param var_range Variance of the range measurements (m^2)
   * @param var_bearing Variance of the bearing measurements (rad^2)
   */
  Measurements(
      std::vector<utias::mrclam::Robot::Measurement> const &measurements,
      double var_range, double var_bearing);

  /**
   * Find the measurement set corresponding to a given timestamp.
   * @param time The timestamp of the desired set of messages.
   * @returns a pointer to a set of measurements. If no mesurements correspond
   * to the given timestamp, a nullptr is return.
   */
  MeasSet *const find(double time);

private:
  /**
   * List of measurements taken by the agent.
   */
  std::vector<MeasSet> measurements_;

  measurementCovariance_t cov_{measurementCovariance_t::Zero()};
};
} // namespace CL::sensors
