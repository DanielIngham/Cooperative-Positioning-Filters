/**
 * @file filter.cpp
 * @brief Implementation of the shared functionality across cooperative
 * localisation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "CL/filters/filter.hpp"

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>
#include <UtiasMrclam/utils/Utils.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CL::filter {

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
        utias::mrclam::utils::getMeasurement(&robot, index)};

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

      if (!utias::mrclam::utils::getMeasurement(robot, k))
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
      if (!utias::mrclam::utils::getMeasurement(robot, k))
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

      if (!utias::mrclam::utils::getMeasurement(&robot, k))
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

} // namespace CL::filter
