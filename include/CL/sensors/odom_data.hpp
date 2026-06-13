#pragma once
#include "CL/common/types.hpp"

/**
 * Data structure housing the list of odometry.
 */
namespace CL::sensors {
struct OdomData {
public:
  OdomData() = delete;

  /**
   * Construct an odometry data instance.
   * @param time Timestamp of the odometry reading.
   * @param vel_f forward velocity in m/s
   * @param vel_w angular velocity in rad/s
   * @param cov 2x2 process noise covariance matrix.
   */
  OdomData(double time, double vel_f, double vel_w, processCovariance_t cov);

  /**
   * Get the timestamp of the odometry input.
   * @returns The timestamp in seconds.
   */
  double time() const;

  /**
   * Get the measured odometry.
   * @returns A constant reference input vector containing the forward and
   * angular velocity.
   */
  input_t const &input() const;

  /**
   * Get the process noise error covariance matrix.
   * @returns A reference to the process noise covariance matrix.
   * @note It was made non-constant to allow for adaptive covariance
   * adjustments.
   */
  processCovariance_t &cov();

private:
  /**
   * Time stamp corresponding to the sensor reading.
   */
  double time_;

  /**
   * Vector containing the forward and angular velocity components of the
   * input.
   */
  input_t input_{input_t::Zero()};

  /**
   * Process noise error covariance matrix.
   */
  processCovariance_t cov_{};
};
} // namespace CL::sensors
