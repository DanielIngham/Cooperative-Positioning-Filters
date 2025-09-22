/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "ekf.h"
#include "filter.h"

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
                     EstimationParameters &estimation_parameters) {

  const double sample_period = data_.getSamplePeriod();

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, estimation_parameters, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(estimation_parameters, sample_period);

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, estimation_parameters, sample_period);

  /* Propagate the estimation error covariance: 3x3 matrix. */
  estimation_parameters.error_covariance =
      estimation_parameters.motion_jacobian *
          estimation_parameters.error_covariance *
          estimation_parameters.motion_jacobian.transpose() +
      estimation_parameters.process_jacobian *
          estimation_parameters.process_noise *
          estimation_parameters.process_jacobian.transpose();
}

#if DECOUPLED
void EKF::correction(EstimationParameters &ego_robot,
                     const EstimationParameters &other_agent) {

  /* Calculate coupled Measurement Jacobian.
   * NOTE: This is only done to reuse functionality between coupled and
   * decoupled approaches.  */
  calculateMeasurementJacobian(ego_robot, other_agent);

  /* Extract the measurment Jacobian in terms of the agent.*/
  Eigen::Matrix<double, 2, 3> agent_measurement_Jacobian =
      ego_robot.measurement_jacobian.topRightCorner<2, 3>();

  /* Calculate joint sensor measurment noise, which is the sum of the measurment
   * noise and the estimate error covariance of the measured agent. */
  measurementCovariance_t joint_sensor_noise =
      ego_robot.measurement_noise + agent_measurement_Jacobian *
                                        other_agent.error_covariance *
                                        agent_measurement_Jacobian.transpose();

  /* Extract the measurement Jacobian in terms of the ego vehicle. */
  Eigen::Matrix<double, 2, 3> ego_measurement_Jacobian =
      ego_robot.measurement_jacobian.topLeftCorner<2, 3>();

  /* Calculate innovation Covariance.  */
  ego_robot.innovation_covariance = ego_measurement_Jacobian *
                                        ego_robot.error_covariance *
                                        ego_measurement_Jacobian.transpose() +
                                    joint_sensor_noise;

  /* Calculate Kalman Gain. */
  Eigen::Matrix<double, 3, 2> kalman_gain =
      ego_robot.error_covariance * ego_measurement_Jacobian.transpose() *
      ego_robot.innovation_covariance.inverse();

  /* Calculate the innovation. */
  measurement_t predicted_measurment = measurementModel(ego_robot, other_agent);

  ego_robot.innovation = ego_robot.measurement - predicted_measurment;
  Data::Robot::normaliseAngle(ego_robot.innovation(BEARING));

  /* Update the state estimate. */
  ego_robot.state_estimate += kalman_gain * ego_robot.innovation;

  /* Update the estimation error covariance.  */
  ego_robot.error_covariance -=
      kalman_gain * ego_robot.innovation_covariance * kalman_gain.transpose();
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
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_agent);

  /* Calculate the innovation: the difference between the measurement
   * and the predicted measurement based on the estimated states of both robots.
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
 * function to increase estimation error covariance of measurements that seem to
 * be outliers.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 */
#if ROBUST
void EKF::robustCorrection(EstimationParameters &ego_robot,
                           const EstimationParameters &other_agent) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t intial_state_estimate = createAugmentedVector(
      ego_robot.state_estimate, other_agent.state_estimate);

  augmentedState_t iterative_state_estimate = intial_state_estimate;

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance = createAugmentedMatrix(
      ego_robot.error_covariance, other_agent.error_covariance);

  /* Calculate the Cholesky Decomposition of the estimation error covariance */
  Eigen::LLT<augmentedCovariance_t> error_cholesky(error_covariance);

  if (error_cholesky.info() != Eigen::Success) {
    std::cout << error_covariance << std::endl;

    throw std::runtime_error("[1] An error has occurred with calculating the "
                             "Cholesky decomposition of "
                             "the estimation error covariance");
  }

  augmentedCovariance_t error_cholesky_matrix = error_cholesky.matrixL();

  /* Calculate the Cholesky Decomposition of the sensor error covariance */
  Eigen::LLT<measurementCovariance_t> measurement_cholesky(
      ego_robot.measurement_noise);

  if (measurement_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the measurement error covariance");
  }

  measurementCovariance_t measurement_cholesky_matrix =
      measurement_cholesky.matrixL();

  /* Calculate measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_agent);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_agent);

  /* Calculate the measurement residual: the difference between the
   * measurement and the calculate measurement based on the estimated states
   * of both robots. */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the bearing residual */
  Data::Robot::normaliseAngle(ego_robot.innovation(BEARING));

  /* Calculate the new robust sensor error covariance. */
  measurementCovariance_t reweighted_measurement_covariance =
      measurement_cholesky_matrix *
      HuberMeasurement(ego_robot.innovation, measurement_thresholds).inverse() *
      measurement_cholesky_matrix.transpose();

  /* Calculate Covariance Innovation. */
  ego_robot.innovation_covariance =
      ego_robot.measurement_jacobian * error_covariance *
          ego_robot.measurement_jacobian.transpose() +
      reweighted_measurement_covariance;

  /* Calculate Kalman Gain. */
  ego_robot.kalman_gain = error_covariance *
                          ego_robot.measurement_jacobian.transpose() *
                          ego_robot.innovation_covariance.inverse();

  /* Update the inititial state estimate. */
  iterative_state_estimate =
      intial_state_estimate + ego_robot.kalman_gain * ego_robot.innovation;

  ego_robot.estimation_residual =
      intial_state_estimate - iterative_state_estimate;

  Data::Robot::normaliseAngle(ego_robot.estimation_residual(ORIENTATION));

  /* Resize matrices back to normal */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Calculate the reweighted error covariance. */
  error_covariance =
      (augmentedCovariance_t::Identity() -
       ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
      error_cholesky_matrix *
      HuberState(ego_robot.estimation_residual, state_thresholds).inverse() *
      error_cholesky_matrix.transpose();

  /* Marginalise the 5x5 back to a 3x3. */
  ego_robot.error_covariance =
      error_covariance.topLeftCorner<total_states, total_states>();
}
#endif // ROBUST

} // namespace Filters
