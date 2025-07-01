/**
 * @file information_filter.cpp
 * @brief Implementation of the Extended Information Filter (information form of
 * the Extended Kalman Filter).
 * @author Daniel Ingham
 * @date 2025-06-21
 */

#include "information_filter.h"
#include "filter.h"

#include <iostream>

/**
 * @brief InformationFilter class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Information filtering.
 * @param[in] data Class containing all robot and landmark data.
 */
InformationFilter::InformationFilter(DataHandler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
InformationFilter::~InformationFilter() {}
/**
 * @brief performs the prediction step of the Information filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] ego_robot The parameters required by the
 * Information filter to perform the prediction step.
 */
#if 0
void InformationFilter::prediction(const Robot::Odometry &odometry,
                                   EstimationParameters &ego_robot) {
  const double sample_period = data_.getSamplePeriod();

  /* Hold onto the previous estimate for use populating the information vector.
   */
  state_t prior_estimate = ego_robot.state_estimate;

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, ego_robot, sample_period);

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, ego_robot, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(ego_robot, sample_period);

  state_t d =
      ego_robot.state_estimate - ego_robot.motion_jacobian * prior_estimate;

  /* Calculate the linearised process covariance and precision. */
  covariance_t linearised_process_covariance =
      ego_robot.process_jacobian * ego_robot.process_noise *
      ego_robot.process_jacobian.transpose();

  // precision_t linearised_process_precision =
  //     linearised_process_covariance.ldlt().solve(precision_t::Identity());
  precision_t linearised_process_precision =
      computePseudoInverse(linearised_process_covariance);

  /* Populate augmented precision. */
  augmentedPrecision_t augmented_precision;

  augmented_precision.topLeftCorner<total_states, total_states>() =
      linearised_process_precision;

  augmented_precision.topRightCorner<total_states, total_states>() =
      -linearised_process_precision * ego_robot.motion_jacobian;

  augmented_precision.bottomLeftCorner<total_states, total_states>() =
      -ego_robot.motion_jacobian.transpose() * linearised_process_precision;

  augmented_precision.bottomRightCorner<total_states, total_states>() =
      ego_robot.motion_jacobian.transpose() * linearised_process_precision *
      ego_robot.motion_jacobian;

  /* Add the contribution of the of the prior precision. */
  augmented_precision.bottomRightCorner<total_states, total_states>() +=
      ego_robot.precision_matrix;

  /* Populate information vector. */
  augmentedInformation_t augmented_information;

  augmented_information.head<total_states>() =
      -linearised_process_precision * d;

  augmented_information.tail<total_states>() =
      ego_robot.motion_jacobian.transpose() * linearised_process_precision * d;

  /* Add the contribution of the prior information. */
  augmented_information.tail<total_states>() += ego_robot.information_vector;

  ego_robot.precision_matrix = marginalise(augmented_precision);

  ego_robot.error_covariance = computePseudoInverse(ego_robot.precision_matrix);

  ego_robot.information_vector =
      marginalise(augmented_information, augmented_precision);

  static bool first = true;

  if (!ego_robot.information_vector.allFinite() && first) {
    first = false;
    std::cout << "Critical Error" << std::endl;
  }

  ego_robot.state_estimate =
      ego_robot.error_covariance * ego_robot.information_vector;
}
#endif
#if 1
void InformationFilter::prediction(const Robot::Odometry &odometry,
                                   EstimationParameters &ego_robot) {

  const double sample_period = data_.getSamplePeriod();

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, ego_robot, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(ego_robot, sample_period);

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, ego_robot, sample_period);

  /* Propagate the estimation information: 3x3 matrix. */
  ego_robot.precision_matrix =
      (ego_robot.motion_jacobian * ego_robot.precision_matrix.inverse() *
           ego_robot.motion_jacobian.transpose() +
       ego_robot.process_jacobian * ego_robot.process_noise *
           ego_robot.process_jacobian.transpose())
          .inverse();

  ego_robot.information_vector =
      ego_robot.precision_matrix * ego_robot.state_estimate;
}
#endif // 0

#ifdef DECOUPLED
void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_agent) {

  /* Calculate the measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_agent);

  Eigen::Matrix<double, 2, 3> ego_measurment_Jacobian =
      ego_robot.measurement_jacobian.topLeftCorner<2, 3>();

  Eigen::Matrix<double, 2, 3> agent_measurment_Jacobian =
      ego_robot.measurement_jacobian.topRightCorner<2, 3>();

  /* Calculate the joint measurment noise. */
  measurementCovariance_t joint_measurement_noise =
      ego_robot.measurement_noise + agent_measurment_Jacobian *
                                        other_agent.precision_matrix.inverse() *
                                        agent_measurment_Jacobian.transpose();

  /* Calculate the measurement residual. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_agent);

  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Calculate the precision contribution */
  precision_t precision_matrix_contribution =
      ego_measurment_Jacobian.transpose() * joint_measurement_noise.inverse() *
      ego_measurment_Jacobian;

  /* Calculate the information contribution */
  information_t information_vector_contribution =
      ego_measurment_Jacobian.transpose() * joint_measurement_noise.inverse() *
      (ego_robot.innovation +
       ego_measurment_Jacobian * ego_robot.state_estimate);

  ego_robot.precision_matrix += precision_matrix_contribution;

  ego_robot.information_vector += information_vector_contribution;

  ego_robot.state_estimate =
      ego_robot.precision_matrix.inverse() * ego_robot.information_vector;
}
#endif

/**
 * @brief Performs Information Filter correct step.
 * @param[in,out] ego_robot The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_agent The robot that was measured by the ego robot.
 * @param[in] robust Flag which determines whether the information and precision
 * should be updated using a robust cost function.
 */
#ifdef COUPLED
void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_agent) {

  /* Create the augmented information vector and precision matrix */
  augmentedInformation_t information_vector = createAugmentedVector(
      ego_robot.information_vector, other_agent.information_vector);

  augmentedPrecision_t precision_matrix = createAugmentedMatrix(
      ego_robot.precision_matrix, other_agent.precision_matrix);

  /* Calculate the augmented estimated state of the system.  */
  augmentedState_t estimated_state =
      computePseudoInverse(precision_matrix) * information_vector;

  /* Calculate measurement Jacobian. */
  calculateMeasurementJacobian(ego_robot, other_agent);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_agent);

  /* Calculate the measurement residual. */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Calculate the precision contribution */
  augmentedPrecision_t precision_matrix_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() * ego_robot.measurement_jacobian;

  /* Calculate the information contribution */
  augmentedInformation_t information_vector_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() *
      (ego_robot.innovation + ego_robot.measurement_jacobian * estimated_state);

  /* Add only the contribution of the of the other agent. */
  precision_matrix_contribution
      .bottomRightCorner<total_states, total_states>() +=
      other_agent.precision_matrix;

  /* Add the information contribution. */
  information_vector_contribution.tail<total_states>() +=
      other_agent.information_vector;

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 6x6 matrix to a 3x3 matrix */
  ego_robot.precision_matrix += marginalise(precision_matrix_contribution);

  /* Marginalise the 6x1 augmented information vector to a 3x1 state vector. */
  ego_robot.information_vector += marginalise(information_vector_contribution,
                                              precision_matrix_contribution);

  /* Retrieve the original state estimate. */
  ego_robot.state_estimate = computePseudoInverse(ego_robot.precision_matrix) *
                             ego_robot.information_vector;
}
#endif // COUPLED
