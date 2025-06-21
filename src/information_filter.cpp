#include "information_filter.h"

InformationFilter::InformationFilter(DataHandler &data) : Filter(data) {}

InformationFilter::~InformationFilter() {}

void InformationFilter::prediction(const Robot::Odometry &odometry,
                                   EstimationParameters &ego_robot) {

  const double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, ego_robot, sample_period);

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, ego_robot, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(ego_robot, sample_period);

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

void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_object,
                                   const bool robust) {

  /* WARN: Robust correction not implemented yet. */
  if (robust) {
    return;
  }

  /* Calculate measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_object);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_object);

  /* Calculate the measurement residual: the difference between the measurement
   * and the calculate measurement based on the estimated states of both robots.
   */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t estimated_state =
      createAugmentedState(ego_robot, other_object);

  /* Calculate the Information contribution */
  augmentedPrecision_t precision_matrix_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() * ego_robot.measurement_jacobian;

  augmentedState_t information_vector_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() *
      (ego_robot.innovation + ego_robot.measurement_jacobian * estimated_state);

  /* Create a temporary augmented matrix  containing the information matrix of
   * both objects. */
  augmentedPrecision_t precision_matrix =
      createAugmentedPrecision(ego_robot, other_object);

  /* Create a temporary augmented vector containing the information vector of
   * both objects. */
  augmentedState_t information_vector = precision_matrix * estimated_state;

  /* Add the information contribution. */
  precision_matrix += precision_matrix_contribution;
  information_vector += information_vector_contribution;

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix */
  ego_robot.precision_matrix = marginalise(precision_matrix);

  ego_robot.information_vector =
      marginalise(information_vector, precision_matrix);

  ego_robot.state_estimate =
      ego_robot.precision_matrix.inverse() * ego_robot.information_vector;
}
