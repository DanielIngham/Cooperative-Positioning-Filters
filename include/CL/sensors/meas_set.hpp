/**
 * @file meas_set.hpp
 */
#pragma once

#include "CL/common/types.hpp"
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
  MeasSet &operator=(MeasSet &&) = delete;
  MeasSet &operator=(const MeasSet &) = delete;
  ~MeasSet() = default;

  /**
   * Construct an instance using the data from the Utias datastructure.
   * @param measurement Utias robot measurement data structure.
   */
  MeasSet(utias::mrclam::Robot::Measurement const &measurement,
          measurementCovariance_t const &cov);

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
  auto begin() const { return meas_set_.begin(); }

  /**
   * Get an iterator to the last element in the set.
   * @returns a read-only (constant) iterator that points to the last element
   * in the set.
   */
  auto end() const { return meas_set_.end(); }

  /**
   * Get the size of the measurement set.
   * @return the size of the %set.
   */
  size_t size() const;

private:
  /** Time stamp of the measurement. */
  double time_;

  /** The set of measurements taken for a given timestamp. */
  std::set<MeasData> meas_set_;

  /** Sensor error covariance matrix. */
  measurementCovariance_t const &cov_{};
};
} // namespace CL::sensors
