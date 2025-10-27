/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "ekf.h"
#include "filter.h"
#include "types.h"

#include <cmath>

namespace Filters {

/**
 * @brief EKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Kalman filtering.
 * @param[in] data Class containing all robot and landmark data.
 */
EKF::EKF(Data::Handler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
EKF::~EKF() {}

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 */
void EKF::prediction(const Data::Robot::Odometry &odometry,
                     EstimationParameters &parameters) {

  const double sample_period{data_.getSamplePeriod()};

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, parameters, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(parameters, sample_period);

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, parameters.state_estimate, sample_period);

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

  const measurementJacobian_t ego_measurement_Jacobian{
      egoMeasurementJacobian(ego, agent)};

  const measurementJacobian_t agent_measurement_Jacobian{
      agentMeasurementJacobian(ego, agent)};

  /* Calculate joint sensor measurment noise, which is the sum of the
   * measurment noise and the estimate error covariance of the measured agent.
   */
  measurementCovariance_t joint_sensor_noise{
      ego.measurement_noise + agent_measurement_Jacobian *
                                  agent.error_covariance *
                                  agent_measurement_Jacobian.transpose()};

  /* Calculate innovation Covariance.  */
  ego.innovation_covariance = ego_measurement_Jacobian * ego.error_covariance *
                                  ego_measurement_Jacobian.transpose() +
                              joint_sensor_noise;

  /* Calculate Kalman Gain. */
  kalmanGain_t kalman_gain{ego.error_covariance *
                           ego_measurement_Jacobian.transpose() *
                           ego.innovation_covariance.inverse()};

  /* Calculate the innovation. */
  measurement_t predicted_measurment{
      measurementModel(ego.state_estimate, agent.state_estimate)};

  ego.innovation = ego.measurement - predicted_measurment;
  Data::Robot::normaliseAngle(ego.innovation(BEARING));

  /* Update the state estimate. */
  ego.state_estimate += kalman_gain * ego.innovation;

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
  calculateMeasurementJacobian(ego_robot, other_agent);

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance = createAugmentedMatrix(
      ego_robot.error_covariance, other_agent.error_covariance);

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
  augmentedState_t state_estimate = createAugmentedVector(
      ego_robot.state_estimate, other_agent.state_estimate);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement{
      measurementModel(ego.state_estimate, agent.state_estimate)};

  /* Calculate the innovation: the difference between the measurement
   * and the predicted measurement based on the estimated states of both
   * robots.
   */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  Data::Robot::normaliseAngle(ego_robot.innovation(BEARING));

  /* Update the state using the measurement. */
  state_estimate += ego_robot.kalman_gain * ego_robot.innovation;

  /* Resize matrices back to normal */
  ego_robot.state_estimate = state_estimate.head<total_states>();

  ego_robot.error_covariance =
      error_covariance.topLeftCorner<total_states, total_states>();
}
#endif // COUPLED

/**
 * @brief A robust version of the correction function that uses the Huber cost
 * function to increase estimation error covariance of measurements that seem
 * to be outliers.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 */

} // namespace Filters
