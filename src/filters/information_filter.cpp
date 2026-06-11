/**
 * @file information_filter.cpp
 * @brief Implementation of the Extended Information Filter (information form of
 * the Extended Kalman Filter).
 * @author Daniel Ingham
 * @date 2025-06-21
 */

#include "CL/filters/information_filter.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include "CL/models/measurement.hpp"
#include "CL/models/process.hpp"
#include "CL/utils/matrix_operations.hpp"
#include "CL/utils/utils.hpp"
#include <iostream>

namespace CL::filter {

EstimationParameters
InformationFilter::prediction(const utias::mrclam::Robot::Odometry &odometry,
                              const EstimationParameters &parameters,
                              double sample_period) {

  EstimationParameters predictive_density{parameters};

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  Models::Process model{odometry, parameters.state_estimate, sample_period};
  const motionJacobian_t motion_jacobian{model.motionJacobian()};
  const processJacobian_t process_jacobian{model.processJacobian()};

  predictive_density.state_estimate = model.predictedState();

  predictive_density.error_covariance =
      motion_jacobian * parameters.precision_matrix.inverse() *
          motion_jacobian.transpose() +
      process_jacobian * parameters.process_noise *
          process_jacobian.transpose();

  predictive_density.precision_matrix =
      predictive_density.error_covariance.inverse();

  predictive_density.information_vector =
      predictive_density.precision_matrix * predictive_density.state_estimate;

  return predictive_density;
}

#ifdef DECOUPLED
/* WARN: This filter approach is unstable due to the inversion of the joint
 * sensor noise matrix. The coupled approach is more robust.   */
void InformationFilter::correction(EstimationParameters &ego,
                                   const EstimationParameters &agent) {
  state_t prior_state{ego.state_estimate};

  Models::Measurement model{ego.state_estimate, agent.state_estimate};
  const measurement_t predicted_measurement{model.predictedMeasurement()};
  const measurementJacobian_t ego_meas_Jacobian{model.egoJacobian()};
  const measurementJacobian_t agent_meas_Jacobian{model.agentJacobian()};

  measurementCovariance_t joint_sensor_noise{
      ego.measurement_noise + agent_meas_Jacobian * agent.error_covariance *
                                  agent_meas_Jacobian.transpose()};

  ego.innovation = ego.measurement - predicted_measurement;
  utils::normaliseAngle(ego.innovation(BEARING));
  ego.innovation += ego_meas_Jacobian * prior_state;

  /* Calculate the precision contribution */
  const precision_t precision_matrix_contribution{
      ego_meas_Jacobian.transpose() * joint_sensor_noise.inverse() *
      ego_meas_Jacobian};

  /* Calculate the information contribution */
  const information_t information_vector_contribution{
      ego_meas_Jacobian.transpose() * joint_sensor_noise.inverse() *
      (ego.innovation)};

  const measurementCovariance_t joint_sensor_precision{
      joint_sensor_noise.inverse()};

  ego.precision_matrix += precision_matrix_contribution;
  ego.error_covariance = ego.precision_matrix.inverse();

  ego.information_vector += information_vector_contribution;

  ego.state_estimate = ego.precision_matrix.inverse() * ego.information_vector;
  utils::normaliseAngle(ego.state_estimate(ORIENTATION));
}
#endif

#ifdef COUPLED
void InformationFilter::correction(EstimationParameters &ego,
                                   const EstimationParameters &agent) {

  /* Create the augmented information vector and precision matrix */
  augmentedInformation_t information_vector{
      MatrixOperations::createAugmentedVector(ego.information_vector,
                                              agent.information_vector)};

  augmentedPrecision_t precision_matrix{MatrixOperations::createAugmentedMatrix(
      ego.precision_matrix, agent.precision_matrix)};

  /* Calculate the augmented estimated state of the system.  */
  augmentedState_t estimated_state{precision_matrix.inverse() *
                                   information_vector};

  /* Calculate measurement Jacobian. */
  Models::Measurement model{ego.state_estimate, agent.state_estimate};
  const augmentedMeasurementJacobian_t measurement_Jacobian{
      model.augmentedJacobian()};

  /* Calculate the measurement residual. */
  ego.innovation = ego.measurement - model.predictedMeasurement();

  /* Normalise the angle residual. */
  utils::normaliseAngle(ego.innovation(BEARING));

  /* Calculate the precision contribution */
  augmentedPrecision_t precision_matrix_contribution{
      measurement_Jacobian.transpose() * ego.measurement_noise.inverse() *
      measurement_Jacobian};

  /* Calculate the information contribution */
  augmentedInformation_t information_vector_contribution{
      measurement_Jacobian.transpose() * ego.measurement_noise.inverse() *
      (ego.innovation + measurement_Jacobian * estimated_state)};

  /* Add only the contribution of the other agent. */
  precision_matrix_contribution
      .bottomRightCorner<total_states, total_states>() +=
      agent.precision_matrix;

  /* Add the information contribution. */
  information_vector_contribution.tail<total_states>() +=
      agent.information_vector;

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 6x6 matrix to a 3x3 matrix */
  ego.precision_matrix +=
      MatrixOperations::marginalise(precision_matrix_contribution);

  /* Marginalise the 6x1 augmented information vector to a 3x1 state vector. */
  ego.information_vector += MatrixOperations::marginalise(
      information_vector_contribution, precision_matrix_contribution);

  /* Retrieve the original state estimate. */
  ego.state_estimate = ego.precision_matrix.inverse() * ego.information_vector;
}
#endif // COUPLED

} // namespace CL::filter
