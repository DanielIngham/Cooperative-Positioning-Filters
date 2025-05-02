/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "ekf.h"
#include "DataHandler/Robot.h"
#include <cmath>

/**
 * @brief EKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Kalman filtering.
 * @param[in] data Class containing all robot data.
 */
EKF::EKF(DataHandler &data) : data_(data) {

  std::vector<Robot> &robots = data_.getRobots();

  /* Populate the Estimation parameters for each robot*/
  for (size_t id = 0; id < data_.getNumberOfRobots(); id++) {
    EstimationParameters initial_parameters;

    /* Assume known prior. This is done by setting the first value of the
     * estimated values to the groundtruth. */
    robots[id].synced.states.push_back(
        Robot::State(robots[id].synced.states.front()));

    initial_parameters.state_estimate << robots[id].synced.states.front().x,
        robots[id].synced.states.front().y,
        robots[id].synced.states.front().orientation;

    /* Populate odometry error covariance matrices */
    initial_parameters.process_noise.diagonal().topRows(total_inputs)
        << robots[id].forward_velocity_error.variance,
        robots[id].angular_velocity_error.variance;

    /* Populate measurement error covariance matrices */
    initial_parameters.measurement_noise.diagonal().topRows(total_inputs)
        << robots[id].range_error.variance,
        robots[id].bearing_error.variance;
  }
}

/**
 * @brief Default destructor.
 */
EKF::~EKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void EKF::peformInference() {
  // std::vector<Robot> &robots = this->data_.getRobots();

  for (size_t k = 1; k < data_.getNumberOfSyncedDatapoints(); k++) {
    for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {
    }
  }
}

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in] prior_state The prior state of the system comprising a x,y
 * coordinate pair and robot orientation.
 * @param[out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 */
void EKF::prediction(const Robot::Odometry &odometry,
                     const Robot::State &prior_state,
                     EstimationParameters &estimation_parameters) {
  double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model. */
  estimation_parameters.state_estimate
      << prior_state.x + odometry.forward_velocity * sample_period *
                             std::cos(prior_state.orientation),
      prior_state.y + odometry.forward_velocity * sample_period *
                          std::sin(prior_state.orientation),
      prior_state.orientation + odometry.angular_velocity * sample_period;

  /* Calculate the Motion Jacobian */
  estimation_parameters.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(prior_state.orientation),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(prior_state.orientation),
      0, 0, 1;

  /* Calculate the process noise Jacobian */
  estimation_parameters.process_jacobian
      << sample_period * std::cos(prior_state.orientation),
      0, sample_period * std::sin(prior_state.orientation), 0, 0, sample_period;

  /* Propagate the estimation error covariance. */
  estimation_parameters.error_covariance =
      estimation_parameters.motion_jacobian *
          estimation_parameters.error_covariance *
          estimation_parameters.motion_jacobian.transpose() +
      estimation_parameters.process_jacobian *
          estimation_parameters.process_noise *
          estimation_parameters.process_jacobian.transpose();
}

/**
 * @brief Perform the correct step.
 */
void EKF::correction(const Robot::Measurement &measurement,
                     const Robot::State &prior_state,
                     EstimationParameters &estimation_parameters,
                     const EstimationParameters &other_robot) {

  double x_difference = estimation_parameters.state_estimate(X, 0) -
                        other_robot.state_estimate(X, 0);
  double y_difference = estimation_parameters.state_estimate(Y, 0) -
                        other_robot.state_estimate(Y, 0);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  /* Create a temporary augmented error covariance matrix to house the
   * error covariance of both the ego vehicle and the measured vehicle:
   * error_covarince =  | P1  0 |
   *                    |  0 P2 |
   * where P1 and P2 are the estimation error covariance of the ego robot and
   * the observed robot respectively.
   */
  Eigen::Matrix<double, 2 * total_states, 2 * total_states> error_covariance;
  error_covariance.setZero();

  error_covariance.topLeftCorner(3, 3) = estimation_parameters.error_covariance;
  error_covariance.bottomRightCorner(3, 3) = other_robot.error_covariance;

  /* Calculate measurement Jacobian */
  estimation_parameters.measurment_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;

  /* Measurement noise Jacobian is identity. No need to calculate. */

  /* TODO: Calculate Innovation */

  /* TODO: Calculate Kalman Gain */

  /* TODO: Update estimation error covariance */

  /* TODO: Correct prediction */
}
