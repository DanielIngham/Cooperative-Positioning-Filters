/**
 * @file information_filter.cpp
 * @brief Implementation of the Extended Information Filter (information form of
 * the Extended Kalman Filter).
 * @author Daniel Ingham
 * @date 2025-06-21
 */

#include "CL/filters/information_filter.hpp"
#include "CL/common/types.hpp"
#include "CL/models/measurement.hpp"
#include "CL/models/process.hpp"
#include "CL/utils/matrix_operations.hpp"

namespace CL::filter {

void InformationFilter::prediction(const Data::Robot::Odometry &odometry,
                                   EstimationParameters &parameters,
                                   double sample_period) {

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  Models::Process model{odometry, parameters.state_estimate, sample_period};
  parameters.state_estimate = model.predictedState();
  const motionJacobian_t &motion_jacobian{model.motionJacobian()};
  const processJacobian_t &process_jacobian{model.processJacobian()};

  /* Propagate the estimation information: 3x3 matrix. */
  parameters.error_covariance = motion_jacobian *
                                    parameters.precision_matrix.inverse() *
                                    motion_jacobian.transpose() +
                                process_jacobian * parameters.process_noise *
                                    process_jacobian.transpose();

  parameters.precision_matrix = parameters.error_covariance.inverse();

  parameters.information_vector =
      parameters.precision_matrix * parameters.state_estimate;
}

#ifdef DECOUPLED
void InformationFilter::correction(EstimationParameters &ego,
                                   const EstimationParameters &agent) {

  const measurementJacobian_t ego_measurement_Jacobian{
      Models::Measurement::egoMeasurementJacobian(ego, agent)};

  const measurementJacobian_t agent_measurement_Jacobian{
      Models::Measurement::agentMeasurementJacobian(ego, agent)};

  /* Calculate the joint measurment noise. */
  const measurementCovariance_t joint_measurement_noise{
      ego.measurement_noise + agent_measurement_Jacobian *
                                  agent.precision_matrix.inverse() *
                                  agent_measurement_Jacobian.transpose()};

  /* Calculate the measurement residual. */
  const measurement_t predicted_measurement{
      Models::Measurement::measurementModel(ego.state_estimate,
                                            agent.state_estimate)};

  ego.innovation = (ego.measurement - predicted_measurement);

  /* Calculate the precision contribution */
  const precision_t precision_matrix_contribution{
      ego_measurement_Jacobian.transpose() * joint_measurement_noise.inverse() *
      ego_measurement_Jacobian};

  /* Calculate the information contribution */
  const information_t information_vector_contribution{
      ego_measurement_Jacobian.transpose() * joint_measurement_noise.inverse() *
      (ego.innovation + ego_measurement_Jacobian * ego.state_estimate)};

  ego.precision_matrix += precision_matrix_contribution;

  ego.information_vector += information_vector_contribution;

  ego.state_estimate = ego.precision_matrix.inverse() * ego.information_vector;
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
  augmentedState_t estimated_state{
      MatrixOperations::computePseudoInverse(precision_matrix) *
      information_vector};

  /* Calculate measurement Jacobian. */
  Models::Measurement::calculateMeasurementJacobian(ego, agent);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement{Models::Measurement::measurementModel(
      ego.state_estimate, agent.state_estimate)};

  /* Calculate the measurement residual. */
  ego.innovation = (ego.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  Data::Robot::normaliseAngle(ego.innovation(BEARING));

  /* Calculate the precision contribution */
  augmentedPrecision_t precision_matrix_contribution{
      ego.measurement_jacobian.transpose() * ego.measurement_noise.inverse() *
      ego.measurement_jacobian};

  /* Calculate the information contribution */
  augmentedInformation_t information_vector_contribution{
      ego.measurement_jacobian.transpose() * ego.measurement_noise.inverse() *
      (ego.innovation + ego.measurement_jacobian * estimated_state)};

  /* Add only the contribution of the of the other agent. */
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
  ego.state_estimate =
      MatrixOperations::computePseudoInverse(ego.precision_matrix) *
      ego.information_vector;
}
#endif // COUPLED

} // namespace CL::filter
