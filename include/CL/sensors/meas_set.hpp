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

  /**
   * Construct an instance using the data from the Utias datastructure.
   * @param measurement Utias robot measurement data structure.
   */
  MeasSet(utias::mrclam::Robot::Measurement const &measurement);

  /**
   * Get message set timestamp.
   * @returns the timestamp of the message.
   */
  double time() const;

  /**
   * Get an iterator to the first element in the set.
   * @returns a read-only (constant) iterator that points to the first element
   * in the set.
   */
  auto begin() const;

  /**
   * Get an iterator to the last element in the set.
   * @returns a read-only (constant) iterator that points to the last element
   * in the set.
   */
  auto end() const;

private:
  /** Time stamp of the measurement. */
  double time_;
  /** The set of measurements taken for a given timestamp. */
  std::set<MeasData> meas_set_;
};
} // namespace CL::sensors
