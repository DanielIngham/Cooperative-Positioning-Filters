/**
 * @file filter.cpp
 * @brief Implementation of the shared functionality across cooperative
 * localisation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#include "filter.h"

#include <DataHandler.h>
#include <chrono>
#include <iostream>

namespace Filter {

/**
 * @brief Assigns fields data based on datahandler input.
 * @param[in] data Class containing all robot data.
 */
Filter::Filter(Data::Handler &data) : data_(data) {

  std::vector<Data::Robot> &robots = data_.getRobots();

  /* Populate the Estimation parameters for each robot. */
  for (auto &robot : robots) {

    EstimationParameters initial_parameters;
    initial_parameters.id = robot.id;
    initial_parameters.barcode = robot.barcode;

    /* Assume known prior. This is done by setting the first value of the
     * estimated values to the groundtruth. */
    robot.synced.states.front() = robot.groundtruth.states.front();

    /* Initial state: 3x1 Matrix. */
    initial_parameters.state_estimate << robot.synced.states.front().x,
        robot.synced.states.front().y, robot.synced.states.front().orientation;

    /* Populate odometry error covariance matrix: 2x2 matrix. */
    initial_parameters.process_noise.diagonal().topRows(total_inputs)
        << robot.forward_velocity_error.variance,
        robot.angular_velocity_error.variance;

    /* Populate measurement error covariance matrix: 2x2 matrix. */
    initial_parameters.measurement_noise.diagonal().topRows(total_measurements)
        << robot.range_error.variance,
        robot.bearing_error.variance;

    /* Populate the Initial information vector */
    initial_parameters.information_vector =
        initial_parameters.precision_matrix * initial_parameters.state_estimate;

    robot_parameters.push_back(initial_parameters);
  }

  /* Populate the estimation parameters for each landmark. */
  std::vector<Data::Landmark> landmarks = data_.getLandmarks();

  for (const auto &landmark : landmarks) {
    EstimationParameters initial_parameters;

    initial_parameters.id = landmark.id;
    initial_parameters.barcode = landmark.barcode;

    initial_parameters.state_estimate << landmark.x, landmark.y, 0.0;

    /* The landmark only has two states: x and y coordintate.
     * NOTE: Although the landmark only has two states, the same data structure
     * is used for both the robot and landmark for compatibility with the
     * EFK::correction function, since only the x and y coordinate and thier
     * corresponding error covariances are used for the measurement update. */
    initial_parameters.error_covariance.diagonal().topRows(total_states - 1)
        << landmark.x_std_dev * landmark.x_std_dev,
        landmark.y_std_dev * landmark.y_std_dev;

    /* NOTE: The landmarks don't have an estimate for thier orientation. So the
     * third diagonal element in the precsion matrix is never used. It is kept
     * for the generality of the implementation of the measurement update step:
     * both robot and landmarks use the same data structure, so they can be used
     * in the same measurement correction function.
     */
    initial_parameters.precision_matrix =
        initial_parameters.error_covariance.inverse();

    initial_parameters.precision_matrix(ORIENTATION, ORIENTATION) = 1;

    initial_parameters.information_vector =
        initial_parameters.precision_matrix * initial_parameters.state_estimate;

    landmark_parameters.push_back(initial_parameters);
  }
}

Filter::~Filter() = default;

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void Filter::performInference() {

  /* Get the total number of syncede datapoints in the data set extracted. */
  size_t total_datapoints = data_.getNumberOfSyncedDatapoints();

  /* Start the timer for measuring the execution time of a child filter. */
  auto timer_start = std::chrono::high_resolution_clock::now();

  for (size_t k{1}; k < total_datapoints; k++) {
    std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints
              << " %" << std::flush;

    /* Perform prediction for each robot using odometry values. */
    std::vector<Data::Robot> &robots = this->data_.getRobots();

    for (unsigned short id{}; id < data_.getNumberOfRobots(); id++) {

      Data::Robot::Odometry odometry = {
          .time = robots[id].synced.odometry[k].time,
          .forward_velocity = robots[id].synced.odometry[k].forward_velocity,
          .angular_velocity = robots[id].synced.odometry[k].angular_velocity,
      };

      prediction(odometry, robot_parameters[id]);

      double normalised_angle =
          robot_parameters[id].state_estimate(ORIENTATION);

      normaliseAngle(normalised_angle);

      /* Update the robot state data structure. */
      robots[id].synced.states[k] = Data::Robot::State{
          .time = robots[id].groundtruth.states[k].time,
          .x = robot_parameters[id].state_estimate(X),
          .y = robot_parameters[id].state_estimate(Y),
          .orientation = normalised_angle,
      };
    }

#ifdef MEASUREMENT_UPDATE
    /* If a measurements are available, loop through each measurement
     * and update the estimate. */
    for (unsigned short id{}; id < data_.getNumberOfRobots(); id++) {

      /* Loop through the measurements taken and perform the measurement
       * update for each robot.
       * NOTE: This operation uses the assumption that the measurements fo the
       * indpendent robots/landmarks are independent of one another.
       */
      const Data::Robot::Measurement *current_measurement =
          Data::Handler::getMeasurement(&robots[id], k);

      if (current_measurement == nullptr) {
        continue;
      }

      for (unsigned short j{}; j < current_measurement->subjects.size(); j++) {
        /* Find the subject for whom the barcode belongs to. */
        Data::Handler::Subject subject;

        bool subject_found =
            data_.getSubject(current_measurement->subjects[j], subject);

        if (!subject_found) {
          continue;
        }

        /* Populate the measurement matrix required for the correction step.
         * Remove any noise bias from the measurement.
         */
        robot_parameters[id].measurement << current_measurement->ranges[j],
            current_measurement->bearings[j];

        /* The datahandler first assigns the ID to the robots then the
         * landmarks. Therefore if the ID is less than or equal to the number
         * of robots, then it belongs to a robot, otherwise it belong to a
         * landmark. */
        EstimationParameters *measured_agent;

        if (subject.type == Data::Handler::Subject::Type::ROBOT) {
          measured_agent = &robot_parameters[subject.index];

        } else {
          measured_agent = &landmark_parameters[subject.index];
        }

        correction(robot_parameters[id], *measured_agent);

        double normalised_angle =
            robot_parameters[id].state_estimate(ORIENTATION);
        normaliseAngle(normalised_angle);

        /* Update the robot state data structure. */
        robots[id].synced.states[k].x = robot_parameters[id].state_estimate(X);
        robots[id].synced.states[k].y = robot_parameters[id].state_estimate(Y);
        robots[id].synced.states[k].orientation = normalised_angle;
      }
    }
#endif // 0
  }

  auto timer_end = std::chrono::high_resolution_clock::now();
  auto execution_duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(timer_end -
                                                            timer_start);

  std::cout << '\r' << std::flush;
  std::cout << "\033[1;32mInference Complete:\033[0m ["
            << execution_duration.count() << " ms]" << std::endl;

  /* Calculate the inference error. */
  data_.calculateStateError();
}

/**
 * @brief Huber Measurement cost function,
 * @param [in] measurement_residual The difference between the measurement and
 * the predicted measurement based on the states of the robots.
 * @param [in] tau Tunable parameter \f$\tau\f$ that is used to determine if the
 * residual is too large.
 */
Filter::huberMeasurementWeights_t
Filter::HuberMeasurement(const measurement_t &measurement_residual,
                         const huberMeasurementThresholds_t &tau) {

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  huberMeasurementWeights_t weight_matrix =
      huberMeasurementWeights_t::Identity();

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
 * @param [in] error_residual The difference between the measurement and
 * the predicted measurement based on the states of the robots.
 * @param [in] tau Tunable parameter \f$\tau\f$ that is used to determine if the
 * residual is too large.
 */
Filter::huberStateWeights_t
Filter::HuberState(const augmentedState_t &error_residual,
                   const huberStateThresholds_t &tau) {

  /* Reweight matrix that is used to adjust the covariance of outliers. */
  huberStateWeights_t weight_matrix = huberStateWeights_t::Identity();

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
 * @brief The unicycle motion model used to perform motion predictions.
 *
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The motion model used for robot takes the form:
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
 * Filter::EstimationParameters::process_noise).
 */
void Filter::motionModel(const Data::Robot::Odometry &odometry,
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
}

/**
 * @brief Calculates the Jacobian matrix of the unicycle motion model in terms
 * of the systems states (x,y,orientation).
 *
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The formula used for the calculation of the motion model
 * Jacobian takes the form:
 * \f[ F = \begin{bmatrix} 1 & 0 & -\tilde{v}\Delta t \sin(\theta) \\ 0 & 1
 * & \tilde{v} \Delta t \cos(\theta) \\ 0 & 0 & 1 \end{bmatrix}, \f]
 * where \f$\theta\f$ denotes the heading (orientation) of the ego vehicle;
 * and \f$\tilde{v}\f$ denotes the forward velocity. The forward velocity is
 * a random variable with Gaussian distributed noise \f$\mathcal{N}(0,w)\f$
 * , where \f$w\f$ is defined by the covariance matrix
 * Filter::EstimationParameters.measurement_noise. See Filter::motionModel for
 * information on the motion model from which this was derived.
 */
void Filter::calculateMotionJacobian(
    const Data::Robot::Odometry &odometry,
    EstimationParameters &estimation_parameters, const double sample_period) {

  estimation_parameters.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, 1;
}

/**
 * @brief Calculates the Jacobian matrix of the motion model evaluated in terms
 * of the process inputs.
 *
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The formula used for the calculation of the process noise
 * Jacobian takes the form
 * \f[L = \begin{bmatrix}\Delta t \cos(\theta) & 0 \\ \Delta t \sin(\theta)
 * & 0 \\ 0 & \Delta t \end{bmatrix}, \f] where \f$\Delta t\f$ denotes the
 * sample period; and \f$\theta\f$ denotes the heading (orientation) of the
 * ego robot. measurement_noise. See EKF::prediction for information on the
 * motion model from which this was derived.
 */
void Filter::calculateProcessJacobian(
    EstimationParameters &estimation_parameters, const double sample_period) {

  estimation_parameters.process_jacobian
      << sample_period *
             std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0,
      sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, sample_period;
}

/**
 * @brief Uses the non-linear measurement model to predict what the measurement
 * from the system would be given the state estimates of the ego robot and
 * measured agent.
 *
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 *
 * @details The measusurement model for the measurement taken from ego vehicle
 * \f$i\f$ to agent \f$j\f$ used for the correction step takes the form
 * \f[ \begin{bmatrix} r_{ij}^{(t)} \\ \phi_{ij}^{(t)}\end{bmatrix} =
 * \begin{bmatrix}\sqrt{(x_j^{(t)} - x_i^{(t)})^2 + (y_j^{(t)} - y_i^{(t)})^2} +
 * q_r \\ \text{atan2}\left(\frac{y_j^{(t)}-y_i^{(t)}}{x_j^{(t)}-x_i^{(t)}
 * }\right) - \theta_i^{(t)} + q_\phi\end{bmatrix}, \f] where \f$x\f$ and
 * \f$y\f$ denote the robots coordinates; \f$\theta\f$ denotes the ego robots
 * orientation (heading); and \f$q_r\f$ and \f$q_\omega\f$ denote the Gaussian
 * distributed measurement noise (See
 * Filter::EstimationParameters.measurement_noise).
 */
Filter::measurement_t
Filter::measurementModel(EstimationParameters &ego_robot,
                         const EstimationParameters &other_agent) {

  /* Calculate the terms */
  const double x_difference =
      other_agent.state_estimate(X) - ego_robot.state_estimate(X);

  const double y_difference =
      other_agent.state_estimate(Y) - ego_robot.state_estimate(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  /* Prevent division by zero and floating point precision errors. */
  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  /* Calculate the predicted measurement based on the estimated states. */
  measurement_t predicted_measurement =
      (measurement_t() << std::sqrt((x_difference * x_difference) +
                                    (y_difference * y_difference)),
       std::atan2(y_difference, x_difference) -
           ego_robot.state_estimate(ORIENTATION))
          .finished();

  normaliseAngle(predicted_measurement(BEARING));

  return predicted_measurement;
}

/**
 * @brief Jacobian of the measurement model evaluated in terms of the systems
 * states: x,y, and heading.
 *
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 *
 * @details The formula used for the calculation of the Jacobian of the
 * measurement matrix between ego vehicle \f$i\f$ and measured agent
 * \f$j\f$ take the form
 * \f[ H = \begin{bmatrix} \frac{-\Delta x}{d} & \frac{-\Delta y}{d} & 0 &
 * \frac{\Delta x}{d} & \frac{\Delta y}{d} & 0\\ \frac{\Delta y}{d^2} &
 * \frac{-\Delta x}{d^2} & -1 & \frac{-\Delta y}{d^2} & \frac{\Delta x}{d^2}
 * & 0\end{bmatrix} \f] where \f$\Delta x = x_j - x_i\f$; \f$\Delta y = y_j
 * - y_i\f$; and \f$\Delta d = \sqrt{\Delta x^2 + \Delta y^2}\f$.
 */
void Filter::calculateMeasurementJacobian(
    EstimationParameters &ego_robot, const EstimationParameters &other_agent) {

  const double x_difference =
      other_agent.state_estimate(X) - ego_robot.state_estimate(X);

  const double y_difference =
      other_agent.state_estimate(Y) - ego_robot.state_estimate(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  ego_robot.measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;
}

/**
 * @brief Schur complement-based marginalisation that marginalises a 6x6 matrix
 * into a 3x3 matrix.
 * @param[in] matrix_6d A 6x6 matrix.
 * @returns A 3x3 matrix.
 */
Filter::matrix3D_t Filter::marginalise(const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight =
      matrix_6d.bottomRightCorner<total_states, total_states>();
  matrix3D_t bottomRightMatrixInverse = computePseudoInverse(bottomRight);

  matrix3D_t matrix_3d =
      matrix_6d.topLeftCorner<total_states, total_states>() -
      matrix_6d.topRightCorner<total_states, total_states>() *
          bottomRightMatrixInverse *
          matrix_6d.bottomLeftCorner<total_states, total_states>();

  return matrix_3d;
}

/**
 * @brief Schur complement-based marginalisation that marginalises a 6x1 vector
 * and 6x6 matrix into a 3x1 vector.
 * @param[in] vector_6d A 6x1 vector.
 * @param[in] matrix_6d A 6x6 matrix.
 * @returns A 3x1 vector.
 */
Filter::state_t Filter::marginalise(const vector6D_t &vector_6d,
                                    const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight =
      matrix_6d.bottomRightCorner<total_states, total_states>();

  matrix3D_t bottomRightInverse = computePseudoInverse(bottomRight);

  state_t new_vector = vector_6d.head<total_states>() -
                       matrix_6d.topRightCorner<total_states, total_states>() *
                           bottomRightInverse * vector_6d.tail<total_states>();

  return new_vector;
}

/**
 * @brief Combines the 3 states of the ego vehicle (x,y,orientation), with the
 * state of the agent measured to create a augmented state vector.
 * @param[in] ego_robot the structure containing the estimation parameters of
 * the ego vehicle.
 * @param[in] other_agent the structure containing the estimation parameters of
 * the agent measured by the ego vehicle.
 * @returns A 5x1 state vector.
 */
Filter::augmentedInformation_t
Filter::createAugmentedVector(const state_t &ego_robot,
                              const state_t &other_agent) {

  augmentedInformation_t augmented_vector = augmentedState_t::Zero();

  augmented_vector.head<total_states>() = ego_robot;
  augmented_vector.tail<total_states>() = other_agent;

  return augmented_vector;
}

/**
 * @brief Combines the covariance/precision matrix of the 3 states of the ego
 * vehicle (x,y,orientation), with the covariance/precision of the of the agent
 * measured.
 *
 * @param[in] ego_robot the structure containing the estimation parameters of
 * the ego vehicle.
 * @param[in] other_agent the structure containing the estimation parameters of
 * the agent measured by the ego vehicle.
 *
 * @details Cooperative Localisation (Positioning) involves robots that share
 * thier state and estimation error covariances / precision when one robot
 * measures the other. As a result, the estimation error covariance/precision
 * needs to be augmented from a 3x3 to a 6x6 matrix to house the error
 * covariance of both the ego vehicle (\f$i\f$) and the measured agent
 * (\f$j\f$):
 * \f[\mathbf{P} = \begin{bmatrix} \mathbf{P}_i & \mathbf{0} \\ \mathbf{0} &
 * \mathbf{P}_j \end{bmatrix}, \f] where \f$\mathbf{P}_i\f$ and
 * \f$\mathbf{P}_j\f$ are the estimation error covariance of the ego robot
 * \f$i\f$ and the observed agent \f$j\f$ respectively.
 */
Filter::augmentedCovariance_t
Filter::createAugmentedMatrix(const covariance_t &ego_robot,
                              const covariance_t &other_agent) {

  augmentedCovariance_t matrix = augmentedCovariance_t::Zero();

  matrix.topLeftCorner<total_states, total_states>() = ego_robot;

  matrix.bottomRightCorner<total_states, total_states>() = other_agent;

  return matrix;
}

/**
 * @brief Normalise an angle between \f$(-\pi, \pi]\f$.
 * @param[inout] angle Angle in radians.
 */
void Filter::normaliseAngle(double &angle) {
  angle -= 2.0 * M_PI * floor((angle + M_PI) / (2.0 * M_PI));
}

/**
 * @brief Calculates the normalised residual of the value produced by the sensor
 * measurement and the prior estimate passed through the non-linear measurement
 * model.
 * @param[in] filter The estimation parameters of the filter.
 * @returns The normalised innovation.
 */
Filter::measurement_t
Filter::calculateNormalisedInnovation(const EstimationParameters &filter) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<measurementCovariance_t> innovatation_cholesky(
      filter.innovation_covariance);

  if (innovatation_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  measurementCovariance_t innovatation_cholesky_matrix =
      innovatation_cholesky.matrixL();

  measurement_t normalised_measurement_residual =
      innovatation_cholesky_matrix.inverse() * filter.innovation;

  return normalised_measurement_residual;
}

/**
 * @brief Calculates the normalised residual of the initial estimate and updated
 * estimate.
 * @param[in] filter The estimation parameters of the filter.
 * @returns The normalised estimation residual.
 */
Filter::augmentedState_t Filter::calculateNormalisedEstimationResidual(
    const EstimationParameters &filter) {

  /* Calculate the mean of the estimation residual (innovation). */
  Eigen::LLT<augmentedCovariance_t> error_covariance_cholesky(
      filter.kalman_gain * filter.innovation_covariance *
          filter.kalman_gain.transpose() +
      augmentedCovariance_t::Identity() * 1e-3);

  if (error_covariance_cholesky.info() != Eigen::Success) {

    throw std::runtime_error("[4] An error has occurred with calculating the "
                             "Cholesky decomposition of "
                             "the estimation error covariance");
  }

  augmentedCovariance_t error_covariance_cholesky_matrix =
      error_covariance_cholesky.matrixL();

  augmentedState_t normalised_error_residual =
      error_covariance_cholesky_matrix.inverse() * filter.estimation_residual;

  return normalised_error_residual;
}

Filter::matrix3D_t Filter::computePseudoInverse(const matrix3D_t &matrix_3d) {

  // Compute pseudo-inverse using SVD
  Eigen::JacobiSVD<matrix3D_t> svd(matrix_3d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV);

  // Set tolerance for singular values (adjust as needed)
  double tolerance = 1e-10;

  // Get singular values and compute pseudo-inverse
  auto singular_values = svd.singularValues();
  Eigen::VectorXd singular_values_inv(singular_values.size());

  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  // Reconstruct pseudo-inverse: V * Σ^+ * U^T
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}

Filter::matrix6D_t Filter::computePseudoInverse(const matrix6D_t &matrix_6d) {

  // Compute pseudo-inverse using SVD
  Eigen::JacobiSVD<matrix6D_t> svd(matrix_6d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV);

  // Set tolerance for singular values (adjust as needed)
  double tolerance = 1e-10;

  // Get singular values and compute pseudo-inverse
  auto singular_values = svd.singularValues();
  Eigen::VectorXd singular_values_inv(singular_values.size());

  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  // Reconstruct pseudo-inverse: V * Σ^+ * U^T
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}
} // namespace Filter
