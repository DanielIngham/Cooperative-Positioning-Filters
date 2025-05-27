/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "iekf.h"
#include "filter.h"
#include <DataHandler/Landmark.h>
#include <DataHandler/Robot.h>
#include <Eigen/Cholesky>
#include <cmath>
#include <iostream>
#include <stdexcept>

/**
 * @brief EKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Kalman filtering.
 * @param[in] data Class containing all robot data.
 */
IEKF::IEKF(DataHandler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
IEKF::~IEKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void IEKF::performInference() {
  std::vector<Robot> &robots = this->data_.getRobots();
  std::vector<Landmark> landmarks = data_.getLandmarks();

  /* Loop through each timestep and perform inference.  */
  std::vector<size_t> measurement_index(data_.getNumberOfRobots(), 1);

  for (size_t k = 1; k < data_.getNumberOfSyncedDatapoints(); k++) {

    /* Perform prediction for each robot using odometry values. */
    for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

      Robot::Odometry odometry(robots[id].synced.odometry[k].time,
                               robots[id].synced.odometry[k].forward_velocity,
                               robots[id].synced.odometry[k].angular_velocity);
      prediction(odometry, robot_parameters[id]);

      /* Update the robot state data structure. */
      robots[id].synced.states.push_back(
          Robot::State(robots[id].groundtruth.states[k].time,
                       robot_parameters[id].state_estimate(X),
                       robot_parameters[id].state_estimate(Y),
                       robot_parameters[id].state_estimate(ORIENTATION)));
    }

    /* If a measurements are available, loop through each measurement
     * and update the estimate. */
    for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

      /* Range check. */
      if (measurement_index[id] >= robots[id].synced.measurements.size()) {
        continue;
      }

      /* Check if a measurement is available. */
      if (std::round(
              (robots[id].synced.measurements[measurement_index[id]].time -
               robots[id].synced.odometry[k].time) *
              10000.0) /
              10000.0 !=
          0.0) {
        continue;
      }

      /* Loop through the measurements taken and perform the measurement
       * update for each robot.
       * NOTE: This operation uses the assumption that the measurements fo the
       * indpendent robots/landmarks are independent of one another.
       */
      const Robot::Measurement &current_measurement =
          robots[id].synced.measurements[measurement_index[id]];

      for (unsigned short j = 0; j < current_measurement.subjects.size(); j++) {
        /* Find the subject for whom the barcode belongs to. */
        int subject_id = data_.getID(current_measurement.subjects[j]);

        if (-1 == subject_id) {
          continue;
        }
        /* Populate the measurement matrix required for the correction step.
         * Remove any noise bias from the measurement.
         */
        robot_parameters[id].measurement << current_measurement.ranges[j],
            current_measurement.bearings[j];

        /* The datahandler first assigns the ID to the robots then the
         * landmarks. Therefore if the ID is less than or equal to the number
         * of robots, then it belongs to a robot, otherwise it belong to a
         * landmark. */
        EstimationParameters measured_object;

        if (subject_id <= data_.getNumberOfRobots()) {
          unsigned short index = subject_id - 1;
          measured_object = robot_parameters[index];

        } else {
          unsigned short index = subject_id - data_.getNumberOfRobots() - 1;
          measured_object = landmark_parameters[index];
        }

        // correction(robot_parameters[id], measured_object);
        robustCorrection(robot_parameters[id], measured_object);

        /* Update the robot state data structure. */
        robots[id].synced.states[k].x = robot_parameters[id].state_estimate(X);
        robots[id].synced.states[k].y = robot_parameters[id].state_estimate(Y);
        robots[id].synced.states[k].orientation =
            robot_parameters[id].state_estimate(ORIENTATION);

        double error = robots[id].synced.states[k].orientation -
                       robots[id].groundtruth.states[k].orientation;

        normaliseAngle(error);

        if (std::abs(error) > 1.5) {
          std::cout << id + 1 << "(" << measured_object.barcode << ")"
                    << ": " << k << "( " << robots[id].synced.states[k].time
                    << ") : " << robots[id].synced.states[k].orientation
                    << " - " << robots[id].groundtruth.states[k].orientation
                    << '\t' << robots[id].synced.states[k - 1].orientation
                    << " - " << robots[id].groundtruth.states[k - 1].orientation
                    << std::endl;

          std::cout << robots[id].synced.odometry[k - 1].forward_velocity << " "
                    << robots[id].synced.odometry[k - 1].angular_velocity
                    << std::endl;
        }
      }
      measurement_index[id] += 1;
    }
  }

  /* Calculate the inference error. */
  for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {
    robots[id].calculateStateError();
  }
}

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 *
 * @details The motion model used for the extended kalman filter prediction take
 * the form
 * \f[\begin{bmatrix} x_i^{(t+1)} \\  y_i^{(t+1)}
 * \\ \theta_i^{(t+1)}\end{bmatrix} = \begin{bmatrix} x_i^{(t)} +
 * \tilde{v}_i^{(t)}\Delta t\cos(\theta_i^{(t)}) \\ y_i^{(t)} +
 * \tilde{v}_i^{(t)}\Delta t\sin(\theta_i^{(t)})\\ \theta_i^{(t)} +
 * \tilde{\omega}_i^{(t)}\Delta t. \end{bmatrix}, \f] where \f$i\f$ denotes the
 * robots ID; \f$t\f$ denotes the current timestep; \f$\Delta t\f$ denotes the
 * sample period; \f$x\f$ and \f$y\f$ are the robots coordinates; \f$\theta\f$
 * denotes the robots heading (orientation); \f$\tilde{v}_i\f$ denotes the
 * forward velocity input; and \f$\tilde{\omega}\f$ denotes the angular velocity
 * input. Both \f$\tilde{v}_i\f$ and \f$\tilde{\omega}_i\f$ are normally
 * distributed random variables \f$\mathcal{N}(0,w)\f$ (see
 * EKF::EstimationParameters::process_noise).
 */
void IEKF::prediction(const Robot::Odometry &odometry,
                      EstimationParameters &estimation_parameters) {

  double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model: 3x1 matrix. */
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

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  estimation_parameters.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, 1;

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  estimation_parameters.process_jacobian
      << sample_period *
             std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0,
      sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, sample_period;

  /* Propagate the estimation error covariance: 3x3 matrix. */
  estimation_parameters.error_covariance =
      estimation_parameters.motion_jacobian *
          estimation_parameters.error_covariance *
          estimation_parameters.motion_jacobian.transpose() +
      estimation_parameters.process_jacobian *
          estimation_parameters.process_noise *
          estimation_parameters.process_jacobian.transpose();
}

/**
 * @brief Performs the Extended Kalman correct step.
 * @param[in,out] ego_robot The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_object The robot that was measured by the ego robot.
 * @details The measusurement model for the measurement taken from ego vehicle
 * \f$i\f$ to vehicle \f$j\f$ used for the correction step takes the form
 * \f[ \begin{bmatrix} r_{ij}^{(t)} \\ \phi_{ij}^{(t)}\end{bmatrix} =
 * \begin{bmatrix}\sqrt{(x_j^{(t)} - x_i^{(t)})^2 + (y_j^{(t)} - y_i^{(t)})^2} +
 * q_r \\ \text{atan2}\left(\frac{y_j^{(t)}-y_i^{(t)}}{x_j^{(t)}-x_i^{(t)}
 * }\right) - \theta_i^{(t)} + q_\phi\end{bmatrix}, \f] where \f$x\f$ and
 * \f$y\f$ denote the robots coordinates; \f$\theta\f$ denotes the ego robots
 * orientation (heading); and \f$q_r\f$ and \f$q_\omega\f$ denote the Gaussian
 * distributed measurement noise (See
 * EKF::EstimationParameters.measurement_noise).
 *
 * @note Cooperative Localisation (Positioning) involves robots that share thier
 * state and estimation error covariances when one robot measures the other. As
 * a result, the estimation error covariance needs to be augmented from a 3x3 to
 * a 6x6 matrix to house the error covariance of both the ego vehicle (\f$i\f$)
 * and the measured vehicle (\f$j\f$):
 * \f[\mathbf{P} = \begin{bmatrix} \mathbf{P}_i & \mathbf{0} \\ \mathbf{0} &
 * \mathbf{P}_j \end{bmatrix}, \f] where \f$\mathbf{P}_i\f$ and
 * \f$\mathbf{P}_j\f$ are the estimation error covariance of the ego robot
 * \f$i\f$ and the observed robot \f$j\f$ respectively.
 */
void IEKF::correction(EstimationParameters &ego_robot,
                      const EstimationParameters &other_object) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  Eigen::Matrix<double, 2 + total_states, 1> intial_state_estimate;
  intial_state_estimate.head<total_states>() = ego_robot.state_estimate;
  intial_state_estimate.tail<total_states - 1>() =
      other_object.state_estimate.head<total_states - 1>();

  Eigen::Matrix<double, 2 + total_states, 1> iterative_state_estimate =
      intial_state_estimate.head<total_states + 2>();

  /* Create and populate new 5x5 error covariance matrix. */
  Eigen::Matrix<double, 2 + total_states, 2 + total_states> error_covariance;
  error_covariance.setZero();

  error_covariance.topLeftCorner<3, 3>() = ego_robot.error_covariance;

  error_covariance.bottomRightCorner<2, 2>() =
      other_object.error_covariance.topLeftCorner<2, 2>();

  /* Perform the iterative update.  */
  const unsigned short max_iterations = 50;

  // std::cout << "New Sensor" << std::endl;
  for (int i = 0; i < max_iterations; i++) {
    /* Calculate measurement Jacobian */
    double x_difference = iterative_state_estimate(X + total_states) -
                          iterative_state_estimate(X);

    double y_difference = iterative_state_estimate(Y + total_states) -
                          iterative_state_estimate(Y);

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
    /* NOTE: Measurement noise Jacobian is identity. No need to calculate. */

    /* Calculate Covariance Innovation: */
    ego_robot.innovation = ego_robot.measurement_jacobian * error_covariance *
                               ego_robot.measurement_jacobian.transpose() +
                           ego_robot.measurement_noise;

    /* Calculate Kalman Gain */
    ego_robot.kalman_gain = error_covariance *
                            ego_robot.measurement_jacobian.transpose() *
                            ego_robot.innovation.inverse();

    /* Populate the predicted measurement matrix. */
    Eigen::Matrix<double, total_measurements, 1> predicted_measurement;

    predicted_measurement << std::sqrt((x_difference * x_difference) +
                                       (y_difference * y_difference)),
        std::atan2(y_difference, x_difference) -
            iterative_state_estimate(ORIENTATION);

    /* Calculate the measurement residual: the difference between the
     * measurement and the calculate measurement based on the estimated states
     * of both robots.
     */
    Eigen::Matrix<double, total_measurements, 1> measurement_residual =
        (ego_robot.measurement - predicted_measurement);

    /* Normalise the bearing residual */
    normaliseAngle(measurement_residual(BEARING));

    // double nis = measurement_residual.transpose() *
    //              ego_robot.innovation.inverse() * measurement_residual;

    // if (nis > 9.21)
    //   std::cout << nis << std::endl;

    Eigen::Matrix<double, 2 + total_states, 1> old_estimate =
        iterative_state_estimate;

    iterative_state_estimate =
        intial_state_estimate +
        ego_robot.kalman_gain *
            (measurement_residual -
             ego_robot.measurement_jacobian *
                 (intial_state_estimate - iterative_state_estimate));

    /* Break if the change between iterations converges */
    double change = (iterative_state_estimate - old_estimate).norm();
    if (change < 1e-8) {
      break;
    }
  }

  /* Resize matrices back to normal */
  /* NOTE: The resize means all information regarding the observed robots
   * states is lost. This operates on the assumptoin that the error covariance
   * is uncorrelated, which may lead to the estimator being over confident in
   * bad estimates. */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(ego_robot.state_estimate(ORIENTATION));

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix by incorporating the
   * marginalising the contributions of the error covariance from the other
   * robot states into the covariance of the ego robot. */
  /* Update estimation error covariance */
  error_covariance -= ego_robot.kalman_gain * ego_robot.innovation *
                      ego_robot.kalman_gain.transpose();

  ego_robot.error_covariance =
      error_covariance.topLeftCorner<total_states, total_states>() -
      error_covariance.topRightCorner<total_states, total_states - 1>() *
          error_covariance
              .bottomRightCorner<total_states - 1, total_states - 1>()
              .inverse() *
          error_covariance.bottomLeftCorner<total_states - 1, total_states>();
}

void IEKF::robustCorrection(EstimationParameters &ego_robot,
                            const EstimationParameters &other_object) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  Eigen::Matrix<double, 2 + total_states, 1> intial_state_estimate;
  intial_state_estimate.head<total_states>() = ego_robot.state_estimate;
  intial_state_estimate.tail<total_states - 1>() =
      other_object.state_estimate.head<total_states - 1>();

  Eigen::Matrix<double, 2 + total_states, 1> iterative_state_estimate =
      intial_state_estimate.head<total_states + 2>();

  /* Create and populate new 5x5 error covariance matrix. */
  Eigen::Matrix<double, 2 + total_states, 2 + total_states> error_covariance;
  error_covariance.setZero();

  error_covariance.topLeftCorner<3, 3>() = ego_robot.error_covariance;

  error_covariance.bottomRightCorner<2, 2>() =
      other_object.error_covariance.topLeftCorner<2, 2>();

  /* TODO: Calculate the Cholesky Decomposition of the estimation error
   * covariance */
  Eigen::LLT<Eigen::Matrix<double, 5, 5>> error_cholesky(error_covariance);

  if (error_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the estimation error covariance");
  }
  Eigen::Matrix<double, 5, 5> error_cholesky_matrix = error_cholesky.matrixL();

  /* TODO: Calculate the Cholesky Decomposition of the sensor error
   * covariance */
  Eigen::LLT<Eigen::Matrix<double, 2, 2>> measurement_cholesky(
      ego_robot.measurement_noise);

  if (measurement_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the measurement error covariance");
  }

  Eigen::Matrix<double, 2, 2> measurement_cholesky_matrix =
      measurement_cholesky.matrixL();

  /* Perform the iterative update.  */
  const unsigned short max_iterations = 50;

  for (int i = 0; i < max_iterations; i++) {
    /* Calculate measurement Jacobian */
    double x_difference = iterative_state_estimate(X + total_states) -
                          iterative_state_estimate(X);

    double y_difference = iterative_state_estimate(Y + total_states) -
                          iterative_state_estimate(Y);

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
    /* NOTE: Measurement noise Jacobian is identity. No need to calculate. */

    /* Populate the predicted measurement matrix. */
    Eigen::Matrix<double, total_measurements, 1> predicted_measurement;
    predicted_measurement << std::sqrt((x_difference * x_difference) +
                                       (y_difference * y_difference)),
        std::atan2(y_difference, x_difference) -
            iterative_state_estimate(ORIENTATION);

    /* Calculate the measurement residual: the difference between the
     * measurement and the calculate measurement based on the estimated states
     * of both robots.
     */
    Eigen::Matrix<double, total_measurements, 1> measurement_residual =
        (ego_robot.measurement - predicted_measurement);

    /* Normalise the bearing residual */
    normaliseAngle(measurement_residual(BEARING));

    /* TODO: Calculate the new robust estimation error covariance. */
    Eigen::Matrix<double, 2 + total_states, 2 + total_states>
        reweighted_error_covariance =
            error_cholesky_matrix *
            HuberState(intial_state_estimate - iterative_state_estimate)
                .inverse() *
            error_cholesky_matrix.transpose();

    /* TODO: Calculate the new robust sensor error covariance. */
    Eigen::Matrix<double, total_measurements, total_measurements>
        reweighted_measurement_covariance =
            measurement_cholesky_matrix *
            HuberMeasurement(measurement_residual).inverse() *
            measurement_cholesky_matrix.transpose();

    /* Calculate Covariance Innovation: */
    ego_robot.innovation = ego_robot.measurement_jacobian *
                               reweighted_error_covariance *
                               ego_robot.measurement_jacobian.transpose() +
                           reweighted_measurement_covariance;

    /* Calculate Kalman Gain */
    ego_robot.kalman_gain = reweighted_error_covariance *
                            ego_robot.measurement_jacobian.transpose() *
                            ego_robot.innovation.inverse();

    Eigen::Matrix<double, 2 + total_states, 1> old_estimate =
        iterative_state_estimate;

    iterative_state_estimate =
        intial_state_estimate +
        ego_robot.kalman_gain *
            (measurement_residual -
             ego_robot.measurement_jacobian *
                 (intial_state_estimate - iterative_state_estimate));

    /* Break if the change between iterations converges */
    double change = (iterative_state_estimate - old_estimate).norm();
    if (change < 1e-8) {
      break;
    }
  }

  /* Resize matrices back to normal */
  /* NOTE: The resize means all information regarding the observed robots
   * states is lost. This operates on the assumptoin that the error covariance
   * is uncorrelated, which may lead to the estimator being over confident in
   * bad estimates. */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(ego_robot.state_estimate(ORIENTATION));

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix by incorporating the
   * marginalising the contributions of the error covariance from the other
   * robot states into the covariance of the ego robot. */
  /* Update estimation error covariance */
  error_covariance =
      (Eigen::Matrix<double, 5, 5>::Identity() -
       ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
      error_cholesky_matrix *
      HuberState(intial_state_estimate - iterative_state_estimate).inverse() *
      error_cholesky_matrix.transpose();

  ego_robot.error_covariance =
      error_covariance.topLeftCorner<total_states, total_states>() -
      error_covariance.topRightCorner<total_states, total_states - 1>() *
          error_covariance
              .bottomRightCorner<total_states - 1, total_states - 1>()
              .inverse() *
          error_covariance.bottomLeftCorner<total_states - 1, total_states>();
}

Eigen::Matrix<double, Filter::total_measurements, Filter::total_measurements>
IEKF::HuberMeasurement(
    const Eigen::Matrix<double, total_measurements, 1> &measurement_residual) {
  // std::cout << "TEST" << std::endl;

  /* Tunable parameter tau that is used to determin if the residual is too large
   * to be an inlier. */
  double tau = 0.5;

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  Eigen::Matrix<double, 2, 2> weight_matrix;
  weight_matrix.setIdentity();

  /* Loop through each of the measurements and perform the huber reweighting if
   * the residual is larger than the parameter tau. */
  bool changed = false;
  for (unsigned short i = 0; i < total_measurements; i++) {

    if (std::abs(measurement_residual(i)) >= tau) {
      changed = true;
      weight_matrix(i, i) = std::abs(tau / measurement_residual(i));
    }
  }
  if (changed) {
    // std::cout << weight_matrix << "\n" << std::endl;
  }

  return weight_matrix;
}

Eigen::Matrix<double, Filter::total_states + 2, Filter::total_states + 2>
IEKF::HuberState(
    const Eigen::Matrix<double, total_states + 2, 1> &error_residual) {

  /* Tunable parameter tau that is used to determine if the residual is too
   * large to be an inlier. */
  double tau = 0.5;

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  Eigen::Matrix<double, total_states + 2, total_states + 2> weight_matrix;
  weight_matrix.setIdentity();

  /* Loop through each of the measurements and perform the huber reweighting if
   * the residual is larger than the parameter tau. */
  bool changed = false;
  for (unsigned short i = 0; i < total_states + 2; i++) {

    if (std::abs(error_residual(i)) >= tau) {
      changed = true;
      weight_matrix(i, i) = std::abs(tau / error_residual(i));
    }
  }

  if (changed) {
    // std::cout << "Weight Matrix After" << std::endl;
    // std::cout << weight_matrix << std::endl;
  }

  return weight_matrix;
}

/**
 * @brief Normalise an angle between \f$\pi\f$ and \f$-\pi\f$.
 * @param[inout] angle angle in radians.
 */
void IEKF::normaliseAngle(double &angle) {
  while (angle >= M_PI)
    angle -= 2.0 * M_PI;

  while (angle < -M_PI)
    angle += 2.0 * M_PI;
}
