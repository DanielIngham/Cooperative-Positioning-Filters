/**
 * @file meas_set.hpp
 */
#pragma once

#include "CL/sensors/meas_data.hpp"

#include <UtiasMrclam/agents/Robot.hpp>
#include <set>

namespace CL::sensors {

/**
 * Contains the set of measurements corresponding to a given timestamp.
 */
class MeasSet {
public:
  MeasSet() = delete;
  MeasSet(MeasSet &&) = default;
  MeasSet(const MeasSet &) = default;
  MeasSet &operator=(MeasSet &&) = default;
  MeasSet &operator=(const MeasSet &) = default;
  ~MeasSet() = default;

  MeasSet(utias::mrclam::Robot::Measurement const &measurement);

  double time() const;

  auto begin() const;

  auto end() const;

private:
  double time_;
  std::set<MeasData> meas_set_;
};
} // namespace CL::sensors
