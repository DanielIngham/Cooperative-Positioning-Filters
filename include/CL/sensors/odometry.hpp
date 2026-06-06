/**
 * @file odometry.hpp
 */
#pragma once

#include "CL/common/types.hpp"
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

  struct Data {
    /**
     * Time stamp corresponding to the sensor reading.
     */
    double time_{};

    /**
     * Vector containing the forward and angular velocity components of the
     * input.
     */
    input_t input_{input_t::Zero()};

    /**
     * Process noise error covariance.
     */
    processCovariance_t cov_{};
  };

  /**
   * Creates an instance housing all the parameters that woul
   * @param odometry Odometry reading for sensor.
   * @param forw_vel_var Forward velocity error variance.
   * @param ang_vel_var Angular velocity error variance.
   */
  Odometry(const std::vector<utias::mrclam::Robot::Odometry> &odometry,
           double forw_vel_var, double ang_vel_var);

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
};
} // namespace CL::sensors
