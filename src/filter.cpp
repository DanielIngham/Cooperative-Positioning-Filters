/**
 * @file filter.cpp
 * @brief Implementation of the shared functionality across cooperative
 * localisation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "filter.h"

#include "Agent.h"
#include "estimation_parameters.h"
#include "types.h"

#include <DataHandler.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace Filters {

/**
 * @brief Assigns fields data based on datahandler input.
 * @param[in] data Class containing all robot data.
 */
Filter::Filter(Data::Handler &data) : data_(data) {

  Data::Robot::List &robots{data_.getRobots()};

  /* Populate the Estimation parameters for each robot. */
  for (auto &robot : robots) {
    /* Assume known prior. This is done by setting the first value of the
     * estimated values to the groundtruth. */
    robot.synced.states.front() = robot.groundtruth.states.front();

    const Data::Agent::ID id{robot.id()};

    auto result{robot_parameters.emplace(id, ParameterList{EstimationParameters{
                                                 .id = id,
                                                 .barcode = robot.barcode(),
                                             }})};
    if (!result.second) {
      throw std::runtime_error("Unable to add Robot with ID " + id +
                               " to robot parameters");
    }

    EstimationParameters &parameters{result.first->second.front()};

    /* Initial state: 3x1 Matrix. */
    parameters.state_estimate(X) = robot.synced.states.front().x;
    parameters.state_estimate(Y) = robot.synced.states.front().y;
    parameters.state_estimate(ORIENTATION) =
        robot.synced.states.front().orientation;

    /* Populate odometry error covariance matrix: 2x2 matrix. */
    parameters.process_noise(FORWARD_VELOCITY, FORWARD_VELOCITY) =
        robot.forward_velocity_error.variance;
    parameters.process_noise(ANGULAR_VELOCITY, ANGULAR_VELOCITY) =
        robot.angular_velocity_error.variance;

    /* Populate measurement error covariance matrix: 2x2 matrix. */
    parameters.measurement_noise(RANGE, RANGE) = robot.range_error.variance;
    parameters.measurement_noise(BEARING, BEARING) =
        robot.bearing_error.variance;

    /* Populate the Initial information vector */
    parameters.information_vector =
        parameters.precision_matrix * parameters.state_estimate;
  }

  /* Populate the estimation parameters for each landmark. */
  const Data::Landmark::List &landmarks{data_.getLandmarks()};

  for (const auto &landmark : landmarks) {
    const Data::Agent::ID id{landmark.id()};

    auto result{
        landmark_parameters.emplace(id, EstimationParameters{
                                            .id = id,
                                            .barcode = landmark.barcode(),
                                        })};
    if (!result.second) {
      throw std::runtime_error("Unable to add Landmark with ID " + id +
                               " to robot parameters");
    }

    EstimationParameters &parameters{result.first->second};

    parameters.state_estimate << landmark.x(), landmark.y(), 0.0;

    /* The landmark only has two states: x and y coordintate.
     * NOTE: Although the landmark only has two states, the same data structure
     * is used for both the robot and landmark for compatibility with the
     * EFK::correction function, since only the x and y coordinate and thier
     * corresponding error covariances are used for the measurement update. */
    parameters.error_covariance.diagonal().topRows(total_states - 1)
        << landmark.x_std_dev() * landmark.x_std_dev(),
        landmark.y_std_dev() * landmark.y_std_dev();

    /* NOTE: The landmarks don't have an estimate for thier orientation. So the
     * third diagonal element in the precsion matrix is never used. It is kept
     * for the generality of the implementation of the measurement update step:
     * both robot and landmarks use the same data structure, so they can be used
     * in the same measurement correction function.
     */
    parameters.precision_matrix = parameters.error_covariance.inverse();

    parameters.precision_matrix(ORIENTATION, ORIENTATION) = 1;

    parameters.information_vector =
        parameters.precision_matrix * parameters.state_estimate;
  }
}

Filter::~Filter() = default;

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void Filter::performInference() {

  /* Get the total number of syncede datapoints in the data set extracted. */
  const size_t total_datapoints{data_.getNumberOfSyncedDatapoints()};

  /* Start the timer for measuring the execution time of a child filter. */
  const auto timer_start{std::chrono::high_resolution_clock::now()};

  /* Perform prediction for each robot using odometry values. */
  Data::Robot::List &robots{data_.getRobots()};

  for (size_t k{1}; k < total_datapoints; k++) {
    std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints
              << " %" << std::flush;

    for (auto &robot : robots) {
      ParameterList &parameter_list{robot_parameters.at(robot.id())};
      parameter_list.push_back(parameter_list.back());

      Data::Robot::Odometry &odometry{robot.synced.odometry[k]};
      EstimationParameters &parameters{parameter_list.back()};

      prediction(odometry, parameters);

      double normalised_angle{parameters.state_estimate(ORIENTATION)};
      Data::Robot::normaliseAngle(normalised_angle);

      /* Update the robot state data structure. */
      robot.synced.states[k] = Data::Robot::State{
          .time = robot.groundtruth.states[k].time,
          .x = parameters.state_estimate(X),
          .y = parameters.state_estimate(Y),
          .orientation = normalised_angle,
      };
    }

#ifdef MEASUREMENT_UPDATE
    /* If a measurements are available, loop through each measurement
     * and update the estimate. */
    processMeasurements(robots, k);
#endif // MEASUREMENT_UPDATE
  }

  const auto timer_end{std::chrono::high_resolution_clock::now()};
  const auto execution_duration{
      std::chrono::duration_cast<std::chrono::milliseconds>(timer_end -
                                                            timer_start)};

  std::cout << '\r' << std::flush;
  std::cout << "\033[1;32mInference Complete:\033[0m ["
            << execution_duration.count() << " ms]" << std::endl;

  /* Calculate the inference error. */
  data_.calculateStateError();
}

void Filter::processMeasurements(Data::Robot::List &robots, size_t index) {

  for (auto &robot : robots) {

    /* Loop through the measurements taken and perform the measurement
     * update for each robot.
     * NOTE: This operation uses the assumption that the measurements fo the
     * indpendent robots/landmarks are independent of one another.
     */
    const Data::Robot::Measurement *current_measurement{
        Data::Handler::getMeasurement(&robot, index)};

    if (!current_measurement) {
      continue;
    }

    /* Populate the measurement matrix required for the correction step.
     * Remove any noise bias from the measurement.
     */
    EstimationParameters &parameters{robot_parameters.at(robot.id()).back()};

    for (unsigned short j{}; j < current_measurement->subjects.size(); j++) {

      /* Find the subject for whom the barcode belongs to. */
      const Data::Agent::Barcode &barcode{current_measurement->subjects[j]};
      EstimationParameters const *measured_agent{
          getEstimationParameters(barcode)};

      if (!measured_agent) {
        continue;
      }

      parameters.measurement[RANGE] = current_measurement->ranges.at(j);
      parameters.measurement[BEARING] = current_measurement->bearings.at(j);

      correction(parameters, *measured_agent);

      double normalised_angle{parameters.state_estimate(ORIENTATION)};

      Data::Robot::normaliseAngle(normalised_angle);

      /* Update the robot state data structure. */
      robot.synced.states.at(index) = {
          .time = robot.groundtruth.states[index].time,
          .x = parameters.state_estimate(X),
          .y = parameters.state_estimate(Y),
          .orientation = normalised_angle,
      };
    }
  }
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
measurement_t
Filter::measurementModel(const EstimationParameters &ego_robot,
                         const EstimationParameters &other_agent) {

  /* Calculate the terms */
  const double x_difference{other_agent.state_estimate(X) -
                            ego_robot.state_estimate(X)};

  const double y_difference{other_agent.state_estimate(Y) -
                            ego_robot.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  /* Prevent division by zero and floating point precision errors. */
  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  /* Calculate the predicted measurement based on the estimated states. */
  measurement_t predicted_measurement{
      (measurement_t() << std::sqrt((x_difference * x_difference) +
                                    (y_difference * y_difference)),
       std::atan2(y_difference, x_difference) -
           ego_robot.state_estimate(ORIENTATION))
          .finished()};

  Data::Robot::normaliseAngle(predicted_measurement(BEARING));

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

  const double x_difference{other_agent.state_estimate(X) -
                            ego_robot.state_estimate(X)};

  const double y_difference{other_agent.state_estimate(Y) -
                            ego_robot.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  ego_robot.measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;
}

measurementJacobian_t
Filter::egoMeasurementJacobian(const EstimationParameters &ego,
                               const EstimationParameters &agent) {

  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  measurementJacobian_t jacobian;
  jacobian(0, 0) = -x_difference / denominator;
  jacobian(0, 1) = -y_difference / denominator;
  jacobian(0, 2) = 0;

  jacobian(1, 0) = y_difference / (denominator * denominator);
  jacobian(1, 1) = -x_difference / (denominator * denominator);
  jacobian(1, 2) = -1;
  return jacobian;
}

measurementJacobian_t
Filter::agentMeasurementJacobian(const EstimationParameters &ego,
                                 const EstimationParameters &agent) {
  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  measurementJacobian_t jacobian;
  jacobian(0, 0) = x_difference / denominator;
  jacobian(0, 1) = y_difference / denominator;
  jacobian(0, 2) = 0;

  jacobian(1, 0) = -y_difference / (denominator * denominator);
  jacobian(1, 1) = x_difference / (denominator * denominator);
  jacobian(1, 2) = 0;

  return jacobian;
}

/**
 * @brief Schur complement-based marginalisation that marginalises a 6x6 matrix
 * into a 3x3 matrix.
 * @param[in] matrix_6d A 6x6 matrix.
 * @returns A 3x3 matrix.
 */
matrix3D_t Filter::marginalise(const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight{
      matrix_6d.bottomRightCorner<total_states, total_states>()};
  matrix3D_t bottomRightMatrixInverse{computePseudoInverse(bottomRight)};

  matrix3D_t matrix_3d{
      matrix_6d.topLeftCorner<total_states, total_states>() -
      matrix_6d.topRightCorner<total_states, total_states>() *
          bottomRightMatrixInverse *
          matrix_6d.bottomLeftCorner<total_states, total_states>()};

  return matrix_3d;
}

/**
 * @brief Schur complement-based marginalisation that marginalises a 6x1 vector
 * and 6x6 matrix into a 3x1 vector.
 * @param[in] vector_6d A 6x1 vector.
 * @param[in] matrix_6d A 6x6 matrix.
 * @returns A 3x1 vector.
 */
state_t Filter::marginalise(const vector6D_t &vector_6d,
                            const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight{
      matrix_6d.bottomRightCorner<total_states, total_states>()};

  matrix3D_t bottomRightInverse{computePseudoInverse(bottomRight)};

  state_t new_vector{vector_6d.head<total_states>() -
                     matrix_6d.topRightCorner<total_states, total_states>() *
                         bottomRightInverse * vector_6d.tail<total_states>()};

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
augmentedInformation_t
Filter::createAugmentedVector(const state_t &ego_robot,
                              const state_t &other_agent) {

  augmentedInformation_t augmented_vector{augmentedState_t::Zero()};

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
augmentedCovariance_t
Filter::createAugmentedMatrix(const covariance_t &ego_robot,
                              const covariance_t &other_agent) {

  augmentedCovariance_t matrix{augmentedCovariance_t::Zero()};

  matrix.topLeftCorner<total_states, total_states>() = ego_robot;

  matrix.bottomRightCorner<total_states, total_states>() = other_agent;

  return matrix;
}

/**
 * @brief Calculates the normalised residual of the value produced by the sensor
 * measurement and the prior estimate passed through the non-linear measurement
 * model.
 * @param[in] filter The estimation parameters of the filter.
 * @returns The normalised innovation.
 */
measurement_t
Filter::normaliseInnovation(const measurement_t &innovation,
                            const measurementCovariance_t &covariance) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<measurementCovariance_t> llt{covariance};

  if (llt.info() != Eigen::Success) {
    std::cout << covariance << std::endl;
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  measurementCovariance_t covariance_sqrt{llt.matrixL()};

  measurement_t normalised_innovation{covariance_sqrt.inverse() * innovation};

  return normalised_innovation;
}

measurement_t
Filter::unnormaliseInnovation(const measurement_t &normalised_innovation,
                              const measurementCovariance_t &covariance) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<measurementCovariance_t> llt{covariance};

  if (llt.info() != Eigen::Success) {
    std::cout << covariance << std::endl;
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  measurementCovariance_t covariance_sqrt{llt.matrixL()};

  measurement_t unnormalised_innovation{covariance_sqrt *
                                        normalised_innovation};

  return unnormalised_innovation;
}

/**
 * @brief Calculates the normalised residual of the initial estimate and updated
 * estimate.
 * @param[in] filter The estimation parameters of the filter.
 * @returns The normalised estimation residual.
 */
augmentedState_t Filter::calculateNormalisedEstimationResidual(
    const EstimationParameters &filter) {

  /* Calculate the mean of the estimation residual (innovation). */
  static constexpr double regularisation{1e-3};

  Eigen::LLT<augmentedCovariance_t> error_covariance_cholesky(
      filter.kalman_gain * filter.innovation_covariance *
          filter.kalman_gain.transpose() +
      augmentedCovariance_t::Identity() * regularisation);

  if (error_covariance_cholesky.info() != Eigen::Success) {

    throw std::runtime_error("[4] An error has occurred with calculating the "
                             "Cholesky decomposition of "
                             "the estimation error covariance");
  }

  augmentedCovariance_t error_covariance_cholesky_matrix{
      error_covariance_cholesky.matrixL()};

  augmentedState_t normalised_error_residual{
      error_covariance_cholesky_matrix.inverse() * filter.estimation_residual};

  return normalised_error_residual;
}

matrix3D_t Filter::computePseudoInverse(const matrix3D_t &matrix_3d) {

  /* Compute pseudo-inverse using SVD */
  Eigen::JacobiSVD<matrix3D_t> svd{matrix_3d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV};

  /* Set tolerance for singular values (adjust as needed) */
  static constexpr double tolerance{1e-10};

  /* Get singular values and compute pseudo-inverse */
  auto singular_values{svd.singularValues()};
  Eigen::VectorXd singular_values_inv{singular_values.size()};

  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  /* Reconstruct pseudo-inverse */
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}

matrix6D_t Filter::computePseudoInverse(const matrix6D_t &matrix_6d) {

  /* Compute pseudo-inverse using SVD */
  Eigen::JacobiSVD<matrix6D_t> svd{matrix_6d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV};

  /* Set tolerance for singular values (adjust as needed) */
  static constexpr double tolerance{1e-10};

  /* Get singular values and compute pseudo-inverse */
  auto singular_values{svd.singularValues()};
  Eigen::VectorXd singular_values_inv{singular_values.size()};

  for (int i{}; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  /* Reconstruct pseudo-inverse. */
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}

/**
 *  Finds the latest estimation parameters corresponding to a give barcode.
 *  @param barcode Barcode of the agent of interests.
 *  @returns Pointer to the agents estimation parameters.
 */
EstimationParameters const *
Filter::getEstimationParameters(const Data::Agent::Barcode &barcode) const {

  for (const auto &[id, robot] : robot_parameters) {
    const EstimationParameters &latest_parameters{robot.back()};

    if (latest_parameters.barcode == barcode) {
      return &latest_parameters;
    }
  }

  for (const auto &[id, landmark] : landmark_parameters) {
    if (landmark.barcode == barcode) {
      return &landmark;
    }
  }

  return nullptr;
}

void Filter::writeInnovation() {
  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};
    if (!success)
      throw std::runtime_error("Unable to create directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {
    const std::string filepath{directory + "/Robot_" + id + "_Innovation.dat"};

    /* Open the file in append mode */
    std::ofstream outFile(filepath, std::ios::app);

    if (!outFile)
      throw std::runtime_error("Failed to open file: " + filepath);

    const Data::Robot *robot{&data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {

      if (!data_.getMeasurement(robot, k))
        continue;

      outFile << parameter_list.at(k).innovation[RANGE] << '\t'
              << parameter_list.at(k).innovation[BEARING] << '\n';
    }

    outFile.close();
  }
}

void Filter::writeNormalisedInnovation() {

  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};

    if (!success)
      throw std::runtime_error("Unable to open directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {

    const std::string filepath{directory + "/Robot_" + id +
                               "_Normalised_Innovation.dat"};

    const Data::Robot *robot{&data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {
      if (!data_.getMeasurement(robot, k))
        continue;

      const measurement_t normalised_innovation{
          normaliseInnovation(parameter_list.at(k).innovation,
                              parameter_list.at(k).innovation_covariance)};

      /* Open the file in append mode */
      std::ofstream outFile(filepath, std::ios::app);
      if (!outFile)
        throw std::runtime_error("Failed to open file: " + filepath);

      outFile << normalised_innovation[RANGE] << '\t'
              << normalised_innovation[BEARING] << '\n';

      outFile.close();
    }
  }
}

void Filter::writeNEES() {

  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};
    if (!success)
      throw std::runtime_error("Unable to create directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {
    const std::string file_path{directory + "/Robot_" + id + "_NEES.dat"};

    const Data::Robot &robot{data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {

      if (!data_.getMeasurement(&robot, k))
        continue;

      const EstimationParameters &parameters{parameter_list.at(k)};

      state_t groundtruth;
      groundtruth << robot.groundtruth.states.at(k).x,
          robot.groundtruth.states.at(k).y,
          robot.groundtruth.states.at(k).orientation;

      state_t error{groundtruth - parameters.state_estimate};
      Data::Robot::normaliseAngle(error(ORIENTATION));

      const double nees{error.transpose() *
                        parameters.error_covariance.inverse() * error};

      /* Open the file in append mode */
      std::ofstream outFile(file_path, std::ios::app);

      if (!outFile)
        throw std::runtime_error("Failed to open file " + file_path);

      outFile << nees << '\n';

      outFile.close();
    }
  }
}

} // namespace Filters
