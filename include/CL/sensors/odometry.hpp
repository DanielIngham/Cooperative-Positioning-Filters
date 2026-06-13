/**
 * @file odometry.hpp
 */
#pragma once

#include "CL/common/types.hpp"
#include "CL/sensors/odom_data.hpp"

#include <Eigen/Dense>
#include <UtiasMrclam/agents/Robot.hpp>

namespace CL::sensors {
/**
 * Data structure that houses all the information required by a filter.
 */
class Odometry {
public:
  Odometry() = default;
  Odometry(Odometry &&) = default;
  Odometry(const Odometry &) = default;
  Odometry &operator=(Odometry &&) = default;
  Odometry &operator=(const Odometry &) = default;
  ~Odometry() = default;

  /**
   * Creates an instance housing all the parameters that woul
   * @param odometry Odometry reading for sensor.
   * @param var_fvel Forward velocity error variance.
   * @param var_avel Angular velocity error variance.
   */
  Odometry(const std::vector<utias::mrclam::Robot::Odometry> &odometry,
           double var_fvel, double var_avel);

  /**
   * Get the timestamp of the odometry message at the given index.
   * @param index Sequence number in the list of odometry.
   */
  double timeAt(size_t index) const;

  /**
   * Get the odometry data at a given sequence number.
   * @param index Sequence number.
   * @returns Odometry data structure.
   */
  OdomData const &odomAt(size_t index) const;

private:
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
  processCovariance_t cov_{processCovariance_t::Zero()};

  /**
   * List of odometry data measured from the odometry sensors.
   */
  std::vector<OdomData> data_;
};
} // namespace CL::sensors
