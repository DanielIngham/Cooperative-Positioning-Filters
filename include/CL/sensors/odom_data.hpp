#pragma once
#include "CL/common/types.hpp"

/**
 * Data structure housing the list of odometry.
 */
namespace CL::sensors {

class OdomData {
public:
  OdomData() = delete;
  OdomData(OdomData &&) = default;
  OdomData(const OdomData &) = default;
  OdomData &operator=(OdomData &&) = default;
  OdomData &operator=(const OdomData &) = default;
  ~OdomData() = default;
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
   * Get either a forward or angular velocity value from the odometry.
   * @returns a double containing either the forward or angular velocity
   * odometry measurement taken.
   */
  double const &input(OdomIdx index) const;

  /**
   * Get the measured odometry.
   * @returns A constant reference input vector containing the forward and
   * angular velocity.
   */
  input_t const &input() const;

  /**
   * Get the process noise error covariance matrix.
   * @returns A constant reference to the process noise covariance matrix.
   */
  processCovariance_t const &noiseCov() const;

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
   * @brief Odometry process noise covariance matrix: 2x2 matrix.
   * @details The process noise covariance matrix is defined by the
   * expression:
   * \f[ w = \begin{bmatrix} q_v & 0 \\ 0 & q_\omega \end{bmatrix}, \f] where
   * \f$q_v\f$ denotes the forward velocity noise variance; and
   * \f$q_\omega\f$ denotes the angular velocity noise variance.
   * @note The process noise is assumed to be uncorrelated and therefore the
   * covariance between the forward velocity and the angular velocity is
   * assumed to be zero.
   */
  processCovariance_t cov_{};
};
} // namespace CL::sensors
