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

/**
 * @brief Performs Information Filter correct step.
 * @param[in,out] ego_robot The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_agent The robot that was measured by the ego robot.
 * @param[in] robust Flag which determines whether the information and precision
 * should be updated using a robust cost function.
 */
void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_agent,
                                   const bool robust) {

  /* WARN: Robust correction not implemented yet. */
  if (robust) {
    static bool first_time_called = true;

    if (first_time_called) {

      first_time_called = false;

      std::cerr
          << "Warning: The robust version of the information filter has not "
             "been "
             "implemented yet. The statndard correction step will be applied."
          << std::endl;
    }
  }
  /* Create the augmented information vector and precision matrix */
  augmentedInformation_t information_vector = augmentedInformation_t::Zero();
  information_vector.head<3>() = ego_robot.information_vector;
  information_vector.tail<3>() = other_agent.information_vector;

  augmentedPrecision_t precision_matrix = augmentedPrecision_t::Zero();
  precision_matrix.topLeftCorner<3, 3>() = ego_robot.precision_matrix;
  precision_matrix.bottomRightCorner<3, 3>() = other_agent.precision_matrix;

  /* Calculate the augmented estimated state of the system.  */
  Eigen::Matrix<double, 6, 1> estimated_state =
      precision_matrix.inverse() * information_vector;

  /* Calculate measurement Jacobian. */
  const double x_difference =
      estimated_state(total_states + X) - estimated_state(X);

  const double y_difference =
      estimated_state(total_states + Y) - estimated_state(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  Eigen::Matrix<double, 2, 6> measurement_jacobian;
  measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_agent);

  /* Calculate the measurement residual. */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Calculate the precision contribution */
  Eigen::Matrix<double, 6, 6> precision_matrix_contribution =
      measurement_jacobian.transpose() * ego_robot.measurement_noise.inverse() *
      measurement_jacobian;

  /* Calculate the information contribution */
  Eigen::Matrix<double, 6, 1> information_vector_contribution =
      measurement_jacobian.transpose() * ego_robot.measurement_noise.inverse() *
      (ego_robot.innovation + measurement_jacobian * estimated_state);

  /* Add only the contribution of the of the other agent. */
  precision_matrix_contribution.bottomRightCorner<3, 3>() +=
      other_agent.precision_matrix;

  /* Add the information contribution. */
  information_vector_contribution.tail<3>() +=
      other_agent.information_vector.head<3>();

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 6x6 matrix to a 3x3 matrix */
  ego_robot.precision_matrix += marginalise(precision_matrix_contribution);

  /* Marginalise the 6x1 augmented information vector to a 3x1 state vector. */
  ego_robot.information_vector += marginalise(information_vector_contribution,
                                              precision_matrix_contribution);

  /* Retrieve the original state estimate. */
  ego_robot.state_estimate =
      ego_robot.precision_matrix.inverse() * ego_robot.information_vector;

  normaliseAngle(ego_robot.state_estimate(ORIENTATION));
}
