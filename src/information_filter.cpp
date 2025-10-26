/**
 * @file information_filter.cpp
 * @brief Implementation of the Extended Information Filter (information form of
 * the Extended Kalman Filter).
 * @author Daniel Ingham
 * @date 2025-06-21
 */

#include "information_filter.h"
#include "filter.h"
#include "types.h"

namespace Filters {
/**
 * @brief InformationFilter class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Information filtering.
 * @param[in] data Class containing all robot and landmark data.
 */
InformationFilter::InformationFilter(Data::Handler &data) : Filter(data) {}

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
void InformationFilter::prediction(const Data::Robot::Odometry &odometry,
                                   EstimationParameters &parameters) {

  const double sample_period{data_.getSamplePeriod()};

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, parameters, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(parameters, sample_period);

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, parameters.state_estimate, sample_period);

  /* Propagate the estimation information: 3x3 matrix. */
  parameters.error_covariance =
      parameters.motion_jacobian * parameters.precision_matrix.inverse() *
          parameters.motion_jacobian.transpose() +
      parameters.process_jacobian * parameters.process_noise *
          parameters.process_jacobian.transpose();

  parameters.precision_matrix = parameters.error_covariance.inverse();

  parameters.information_vector =
      parameters.precision_matrix * parameters.state_estimate;
}

#ifdef DECOUPLED
void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_agent) {

  const measurementJacobian_t ego_measurement_Jacobian{
      egoMeasurementJacobian(ego_robot, other_agent)};

  const measurementJacobian_t agent_measurement_Jacobian{
      agentMeasurementJacobian(ego_robot, other_agent)};

  /* Calculate the joint measurment noise. */
  const measurementCovariance_t joint_measurement_noise{
      ego_robot.measurement_noise + agent_measurement_Jacobian *
                                        other_agent.precision_matrix.inverse() *
                                        agent_measurement_Jacobian.transpose()};

  /* Calculate the measurement residual. */
  const measurement_t predicted_measurement{
      measurementModel(ego_robot, other_agent)};

  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Calculate the precision contribution */
  const precision_t precision_matrix_contribution{
      ego_measurement_Jacobian.transpose() * joint_measurement_noise.inverse() *
      ego_measurement_Jacobian};

  /* Calculate the information contribution */
  const information_t information_vector_contribution{
      ego_measurement_Jacobian.transpose() * joint_measurement_noise.inverse() *
      (ego_robot.innovation +
       ego_measurement_Jacobian * ego_robot.state_estimate)};

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
  Data::Robot::normaliseAngle(ego_robot.innovation(BEARING));

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

} // namespace Filters
