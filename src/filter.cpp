#include "filter.h"
#include <DataHandler/DataHandler.h>
#include <iostream>

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
    // initial_parameters.process_noise.diagonal().topRows(total_inputs)
    //     << robots[id].forward_velocity_error.variance,
    //     robots[id].angular_velocity_error.variance;

    initial_parameters.process_noise.diagonal().topRows(total_inputs) << 0.0001,
        0.0001;
    std::cout << " Process Noise " << std::endl;
    std::cout << initial_parameters.process_noise << std::endl;
    std::cout << " " << std::endl;

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
      innovatation_cholesky_matrix.inverse() * filter.measurement_residual;

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
