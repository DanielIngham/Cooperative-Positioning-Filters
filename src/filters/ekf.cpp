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
#include "CL/utils/utils.hpp"

#include <cmath>

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
  Models::Process model{odometry, parameters.state_estimate, sample_period};

  parameters.state_estimate = model.predictedState();
  const motionJacobian_t &motion_jacobian{model.motionJacobian()};
  const processJacobian_t &process_jacobian{model.processJacobian()};

  /* Propagate the estimation error covariance: 3x3 matrix. */
  parameters.error_covariance = motion_jacobian * parameters.error_covariance *
                                    motion_jacobian.transpose() +
                                process_jacobian * parameters.process_noise *
                                    process_jacobian.transpose();
}

#if DECOUPLED
void EKF::correction(EstimationParameters &ego,
                     const EstimationParameters &agent) {

  const measurementJacobian_t ego_measurement_Jacobian{
      Models::Measurement::egoMeasurementJacobian(ego, agent)};

  const measurementJacobian_t agent_measurement_Jacobian{
      Models::Measurement::agentMeasurementJacobian(ego, agent)};

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
  measurement_t predicted_measurment{Models::Measurement::measurementModel(
      ego.state_estimate, agent.state_estimate)};

  ego.innovation = ego.measurement - predicted_measurment;
  utils::normaliseAngle(ego.innovation(BEARING));

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

#if RANGE_ONLY
void EKF::correction(EstimationParameters &ego,
                     const EstimationParameters &agent) {

  const vector3D_t ego_measurement_Jacobian{
      Models::Measurement::egoRangeMeasurementJacobian(ego, agent)};

  const vector3D_t agent_measurement_Jacobian{
      Models::Measurement::agentRangeMeasurementJacobian(ego, agent)};

  /* Calculate joint sensor measurment noise, which is the sum of the
   * measurment noise and the estimate error covariance of the measured agent.
   */
  double ego_measurement_noise{ego.measurement_noise(0)};

  double joint_sensor_noise{ego_measurement_noise +
                            agent_measurement_Jacobian.transpose() *
                                agent.error_covariance *
                                agent_measurement_Jacobian};

  /* Calculate innovation Covariance.  */
  double innovation_covariance{ego_measurement_Jacobian.transpose() *
                                   ego.error_covariance *
                                   ego_measurement_Jacobian +
                               joint_sensor_noise};

  /* Calculate Kalman Gain. */
  vector3D_t kalman_gain{(ego.error_covariance * ego_measurement_Jacobian) /
                         innovation_covariance};

  /* Calculate the innovation. */
  const double predicted_measurement{Models::Measurement::rangeMeasurementModel(
      ego.state_estimate, agent.state_estimate)};

  const double measurement_range{ego.measurement(RANGE)};
  const double innovation{measurement_range - predicted_measurement};

  /* Update the state estimate. */
  ego.state_estimate += kalman_gain * innovation;

  /* Update the estimation error covariance.  */
  ego.error_covariance -=
      kalman_gain * innovation_covariance * kalman_gain.transpose();
}

#endif // BEARING_ONLY

} // namespace CL::filter
