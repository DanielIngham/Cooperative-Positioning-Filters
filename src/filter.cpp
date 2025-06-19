#include "filter.h"
#include <DataHandler/DataHandler.h>

/**
 * @brief Assigns fields data based on datahandler input.
 * @param[in] data Class containing all robot data.
 */
Filter::Filter(DataHandler &data) : data_(data) {

  std::vector<Robot> &robots = data_.getRobots();

  /* Populate the Estimation parameters for each robot. */
  for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

    EstimationParameters initial_parameters;
    initial_parameters.id = robots[id].id;
    initial_parameters.barcode = robots[id].barcode;

    /* Assume known prior. This is done by setting the first value of the
     * estimated values to the groundtruth. */
    robots[id].synced.states.push_back(
        Robot::State(robots[id].groundtruth.states.front()));

    /* Initial state: 3x1 Matrix. */
    initial_parameters.state_estimate << robots[id].synced.states.front().x,
        robots[id].synced.states.front().y,
        robots[id].synced.states.front().orientation;

    /* Populate odometry error covariance matrix: 2x2 matrix. */
    initial_parameters.process_noise.diagonal().topRows(total_inputs)
        << robots[id].forward_velocity_error.variance,
        robots[id].angular_velocity_error.variance;

    /* Populate measurement error covariance matrix: 2x2 matrix. */
    initial_parameters.measurement_noise.diagonal().topRows(total_measurements)
        << robots[id].range_error.variance,
        robots[id].bearing_error.variance;

    /* Populate the Initial information vector */
    initial_parameters.information_vector =
        initial_parameters.precision_matrix * initial_parameters.state_estimate;

    robot_parameters.push_back(initial_parameters);
  }

  /* Populate the estimation parameters for each landmark. */
  std::vector<Landmark> landmarks = data_.getLandmarks();

  for (unsigned short id = 0; id < data_.getNumberOfLandmarks(); id++) {
    EstimationParameters initial_parameters;

    initial_parameters.id = landmarks[id].id;
    initial_parameters.barcode = landmarks[id].barcode;

    initial_parameters.state_estimate << landmarks[id].x, landmarks[id].y, 0.0;

    /* The landmark only has two states: x and y coordintate.
     * NOTE: Although the landmark only has two states, the same data structure
     * is used for both the robot and landmark for compatibility with the
     * EFK::correction function, since only the x and y coordinate and thier
     * corresponding error covariances are used for the measurement update. */
    initial_parameters.error_covariance.diagonal().topRows(total_states - 1)
        << landmarks[id].x_std_dev * landmarks[id].x_std_dev,
        landmarks[id].y_std_dev * landmarks[id].y_std_dev;

    /* NOTE: The landmarks don't have an estimate for thier orientation. So the
     * third diagonal element in the precsion matrix is never used. It is kept
     * for the generality of the implementation of the measurement update step:
     * both robot and landmarks use the same data structure, so they can be used
     * in the same measurement correction function.
     */
    initial_parameters.precision_matrix =
        initial_parameters.error_covariance.inverse();

    landmark_parameters.push_back(initial_parameters);
  }
}

Filter::~Filter() = default;

/**
 * @brief Huber Measurement cost function,
 * @param [in] measurement_residual The difference between the measurement and
 * the predicted measurement based on the states of the robots.
 * @param [in] tau Tunable parameter \f$\tau\f$ that is used to determine if the
 * residual is too large.
 */
Eigen::Matrix<double, Filter::total_measurements, Filter::total_measurements>
Filter::HuberMeasurement(
    const Eigen::Matrix<double, total_measurements, 1> &measurement_residual,
    const Eigen::Matrix<double, total_measurements, 1> &tau) {

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  Eigen::Matrix<double, total_measurements, total_measurements> weight_matrix =
      Eigen::Matrix<double, total_measurements, total_measurements>::Identity();

  /* Loop through each of the measurements and perform the huber reweighting if
   * the residual is larger than the parameter tau. */
  for (unsigned short i = 0; i < total_measurements; i++) {

    if (std::abs(measurement_residual(i)) >= tau(i)) {
      weight_matrix(i, i) = tau(i) / std::abs(measurement_residual(i));
    }
  }

  return weight_matrix;
}

/**
 * @brief Huber State cost function,
 * @param [in] measurement_residual The difference between the measurement and
 * the predicted measurement based on the states of the robots.
 * @param [in] tau Tunable parameter \f$\tau\f$ that is used to determine if the
 * residual is too large.
 */
Eigen::Matrix<double, 2 + Filter::total_states, 2 + Filter::total_states>
Filter::HuberState(
    const Eigen::Matrix<double, 2 + total_states, 1> &error_residual,
    const Eigen::Matrix<double, 2 + total_states, 1> &tau) {

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  Eigen::Matrix<double, 2 + total_states, 2 + total_states> weight_matrix =
      Eigen::Matrix<double, 2 + total_states, 2 + total_states>::Identity();

  /* Loop through each of the measurements and perform the huber reweighting if
   * the residual is larger than the parameter tau. */
  for (unsigned short i = 0; i < 2 + total_states; i++) {

    if (std::abs(error_residual(i)) >= tau(i)) {
      weight_matrix(i, i) = tau(i) / std::abs(error_residual(i));
    }
  }

  return weight_matrix;
}

void Filter::motionModel(const Robot::Odometry &odometry,
                         EstimationParameters &estimation_parameters,
                         const double sample_period) {

  estimation_parameters.state_estimate
      << estimation_parameters.state_estimate(X) +
             odometry.forward_velocity * sample_period *
                 std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      estimation_parameters.state_estimate(Y) +
          odometry.forward_velocity * sample_period *
              std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      estimation_parameters.state_estimate(ORIENTATION) +
          odometry.angular_velocity * sample_period;

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(estimation_parameters.state_estimate(ORIENTATION));
}

void Filter::motionJacobian(const Robot::Odometry &odometry,
                            EstimationParameters &estimation_parameters,
                            const double sample_period) {

  estimation_parameters.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, 1;
}

void Filter::processJacobian(EstimationParameters &estimation_parameters,
                             const double sample_period) {

  estimation_parameters.process_jacobian
      << sample_period *
             std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0,
      sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, sample_period;
}

Eigen::Matrix<double, Filter::total_measurements, 1>
Filter::measurementModel(EstimationParameters &ego_robot,
                         const EstimationParameters &other_object) {

  const double x_difference =
      other_object.state_estimate(X) - ego_robot.state_estimate(X);

  const double y_difference =
      other_object.state_estimate(Y) - ego_robot.state_estimate(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  Eigen::Matrix<double, total_measurements, 1> predicted_measurement;

  predicted_measurement << std::sqrt((x_difference * x_difference) +
                                     (y_difference * y_difference)),
      std::atan2(y_difference, x_difference) -
          ego_robot.state_estimate[ORIENTATION];

  normaliseAngle(predicted_measurement(BEARING));

  return predicted_measurement;
}

void Filter::calculateMeasurementJacobian(
    EstimationParameters &ego_robot, const EstimationParameters &other_object) {

  const double x_difference =
      other_object.state_estimate(X) - ego_robot.state_estimate(X);

  const double y_difference =
      other_object.state_estimate(Y) - ego_robot.state_estimate(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  ego_robot.measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator);
}

/**
 * @brief Schur complement-based error covariance marginalisation.
 * @details This is used to marginalise the 5x5 matrix to a 3x3 matrix by
 * incorporating the marginalising the contributions of the error covariance
 * from the other robot states into the covariance of the ego robot.
 */
Eigen::Matrix<double, Filter::total_states, Filter::total_states>
Filter::marginalise(
    Eigen::Matrix<double, 2 + Filter::total_states, 2 + Filter::total_states>
        matrix_5d) {

  Eigen::Matrix<double, total_states, total_states> matrix_3d =
      matrix_5d.topLeftCorner<total_states, total_states>() -
      matrix_5d.topRightCorner<total_states, total_states - 1>() *
          matrix_5d.bottomRightCorner<total_states - 1, total_states - 1>()
              .inverse() *
          matrix_5d.bottomLeftCorner<total_states - 1, total_states>();

  return matrix_3d;
}

Eigen::Matrix<double, 2 + Filter::total_states, 1>
Filter::createAugmentedState(const EstimationParameters &ego_robot,
                             const EstimationParameters &other_object) {

  Eigen::Matrix<double, 2 + Filter::total_states, 1> state_estimate =
      Eigen::Matrix<double, 2 + Filter::total_states, 1>::Zero();

  state_estimate.head<total_states>() = ego_robot.state_estimate;
  state_estimate.tail<total_states - 1>() =
      other_object.state_estimate.head<total_states - 1>();

  return state_estimate;
}

Eigen::Matrix<double, 2 + Filter::total_states, 2 + Filter::total_states>
Filter::createAugmentedCovariance(const EstimationParameters &ego_robot,
                                  const EstimationParameters &other_object) {
  Eigen::Matrix<double, 2 + Filter::total_states, 2 + Filter::total_states>
      matrix;

  matrix.setZero();

  matrix.topLeftCorner<3, 3>() = ego_robot.error_covariance;

  matrix.bottomRightCorner<2, 2>() =
      other_object.error_covariance.topLeftCorner<2, 2>();

  return matrix;
}

/**
 * @brief Normalise an angle between \f$\pi\f$ and \f$-\pi\f$.
 * @param[inout] angle angle in radians.
 */
void Filter::normaliseAngle(double &angle) {
  while (angle >= M_PI)
    angle -= 2.0 * M_PI;

  while (angle < -M_PI)
    angle += 2.0 * M_PI;
}

Eigen::Matrix<double, Filter::total_measurements, 1>
Filter::calculateNormalisedMeasurementResidual(
    const EstimationParameters &filter) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<Eigen::Matrix<double, total_measurements, total_measurements>>
      innovatation_cholesky(filter.innovation_covariance);

  if (innovatation_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  Eigen::Matrix<double, total_measurements, total_measurements>
      innovatation_cholesky_matrix = innovatation_cholesky.matrixL();

  Eigen::Matrix<double, total_measurements, 1> normalised_measurement_residual =
      innovatation_cholesky_matrix.inverse() * filter.innovation;

  return normalised_measurement_residual;
}

Eigen::Matrix<double, 2 + Filter::total_states, 1>
Filter::calculateNormalisedEstimationResidual(
    const EstimationParameters &filter) {

  /* Calculate the mean of the estimation residual (innovation). */
  Eigen::LLT<Eigen::Matrix<double, 2 + total_states, 2 + total_states>>
      error_covariance_cholesky(filter.kalman_gain *
                                    filter.innovation_covariance *
                                    filter.kalman_gain.transpose() +
                                Eigen::Matrix<double, 2 + total_states,
                                              2 + total_states>::Identity() *
                                    1e-3);

  if (error_covariance_cholesky.info() != Eigen::Success) {

    throw std::runtime_error("[4] An error has occurred with calculating the "
                             "Cholesky decomposition of "
                             "the estimation error covariance");
  }

  Eigen::Matrix<double, 2 + total_states, 2 + total_states>
      error_covariance_cholesky_matrix = error_covariance_cholesky.matrixL();

  Eigen::Matrix<double, 2 + total_states, 1> normalised_error_residual =
      error_covariance_cholesky_matrix.inverse() * filter.estimation_residual;

  return normalised_error_residual;
}
