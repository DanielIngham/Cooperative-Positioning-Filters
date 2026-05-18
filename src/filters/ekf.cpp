/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "CL/filters/ekf.hpp"
#include "CL/common/types.hpp"
#include "CL/models/measurement.hpp"
#include "CL/models/process.hpp"
#include "CL/models/range.hpp"
#include "CL/models/range_bearing.hpp"
#include "CL/utils/utils.hpp"

#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

#ifdef COUPLED
#include "CL/utils/matrix_operations.hpp"
#endif

namespace CL::filter {

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 */
void EKF::prediction(const Data::Robot::Odometry &odometry,
                     EstimationParameters &parameters, double sample_period) {

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  Models::Process::calculateMotionJacobian(odometry, parameters, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  Models::Process::calculateProcessJacobian(parameters, sample_period);

  /* Make the prediction using the motion model: 3x1 matrix. */
  Models::Process::motionModel(odometry, parameters.state_estimate,
                               sample_period);

  /* Propagate the estimation error covariance: 3x3 matrix. */
  parameters.error_covariance =
      parameters.motion_jacobian * parameters.error_covariance *
          parameters.motion_jacobian.transpose() +
      parameters.process_jacobian * parameters.process_noise *
          parameters.process_jacobian.transpose();
}

#if DECOUPLED
void EKF::correction(EstimationParameters &ego,
                     const EstimationParameters &agent) {

  const auto measurement{
      Models::Measurement::generateMeasurement<Models::RangeBearing>(ego,
                                                                     agent)};

  /* The joint sensor measurement noise is the sum of the measurement noise and
   * the estimate error covariance of the measured agent. */
  Eigen::MatrixXd joint_sensor_noise{
      ego.measurement_noise + measurement.getAgentJacobian() *
                                  agent.error_covariance *
                                  measurement.getAgentJacobian().transpose()};

  ego.innovation_covariance = measurement.getEgoJacobian() *
                                  ego.error_covariance *
                                  measurement.getEgoJacobian().transpose() +
                              joint_sensor_noise;

  Eigen::MatrixXd kalman_gain{ego.error_covariance *
                              measurement.getEgoJacobian().transpose() *
                              ego.innovation_covariance.inverse()};

  ego.innovation = ego.measurement - measurement.getPrediction();
  utils::normaliseAngle(ego.innovation(BEARING));

  /* Update the state estimate. */
  ego.state_estimate += kalman_gain * ego.innovation;

  if (ego.state_estimate.hasNaN()) {
    std::cerr << "Ego Jacobian:\n" << measurement.getEgoJacobian() << std::endl;
    std::cerr << "Agent Jacobian:\n"
              << measurement.getAgentJacobian() << std::endl;
    std::cerr << "Innovation Covariance:\n"
              << ego.innovation_covariance << std::endl;
    std::cerr << "Innovation:\n" << ego.innovation << std::endl;
    std::cerr << "Kalman Gain:\n" << ego.kalman_gain << std::endl;
    std::cerr << "Predicted Measurement" << measurement.getPrediction()
              << std::endl;

    throw std::runtime_error("State estimate contains a NaN");
  }

  /* Update the estimation error covariance.  */
  ego.error_covariance -=
      kalman_gain * ego.innovation_covariance * kalman_gain.transpose();
}
#endif // DECOUPLED

/**
 * @brief Performs the Extended Kalman correct step.
 * @param[in,out] ego_robot The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_agent The agent that was measured by the ego robot.
 * @param[in] robust Flag which determines whether the state and covariance
 * should be updated using a robust cost function.
 */
#if COUPLED

void EKF::correction(EstimationParameters &ego_robot,
                     const EstimationParameters &other_agent) {

  /* Calculate measurement Jacobian */
  Models::Measurement::calculateMeasurementJacobian(ego_robot, other_agent);

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      MatrixOperations::createAugmentedMatrix(ego_robot.error_covariance,
                                              other_agent.error_covariance);

  /* Calculate Covariance Innovation: */
  ego_robot.innovation_covariance =
      ego_robot.measurement_jacobian * error_covariance *
          ego_robot.measurement_jacobian.transpose() +
      ego_robot.measurement_noise;

  /* Calculate Kalman Gain */
  ego_robot.kalman_gain = error_covariance *
                          ego_robot.measurement_jacobian.transpose() *
                          ego_robot.innovation_covariance.inverse();

  /* Update estimation error covariance */
  error_covariance -= ego_robot.kalman_gain * ego_robot.innovation_covariance *
                      ego_robot.kalman_gain.transpose();

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t state_estimate = MatrixOperations::createAugmentedVector(
      ego_robot.state_estimate, other_agent.state_estimate);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement{Models::Measurement::measurementModel(
      ego_robot.state_estimate, other_agent.state_estimate)};

  /* Calculate the innovation: the difference between the measurement
   * and the predicted measurement based on the estimated states of both
   * robots.
   */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  utils::normaliseAngle(ego_robot.innovation(BEARING));

  /* Update the state using the measurement. */
  state_estimate += ego_robot.kalman_gain * ego_robot.innovation;

  /* Resize matrices back to normal */
  ego_robot.state_estimate = state_estimate.head<total_states>();

  ego_robot.error_covariance =
      error_covariance.topLeftCorner<total_states, total_states>();
}
#endif // COUPLED

} // namespace CL::filter
