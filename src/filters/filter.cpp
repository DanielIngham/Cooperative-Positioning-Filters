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
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CL::filter {

Filter::~Filter() = default;

void Filter::performInference() {

  /* Start the timer for measuring the execution time of a child filter. */
  const auto timer_start{std::chrono::high_resolution_clock::now()};

  for (size_t k{1}; k < total_datapoints_; k++) {
    std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints_
              << " %" << std::flush;

    for (auto &robot : robots_) {
      ParameterList &parameter_list{robot_parameters.at(robot.id())};
      parameter_list.push_back(parameter_list.back());

      Data::Robot::Odometry &odometry{robot.synced.odometry[k]};
      EstimationParameters &parameters{parameter_list.back()};

      prediction(odometry, parameters);

      double normalised_angle{parameters.state_estimate(ORIENTATION)};
      utias::mrclam::utils::normaliseAngle(normalised_angle);

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
    processMeasurements(robots_, k);
#endif // MEASUREMENT_UPDATE
  }

  const auto timer_end{std::chrono::high_resolution_clock::now()};
  const auto execution_duration{
      std::chrono::duration_cast<std::chrono::milliseconds>(timer_end -
                                                            timer_start)};

  std::cout << '\r' << std::flush;
  std::cout << "\033[1;32mInference Complete:\033[0m ["
            << execution_duration.count() << " ms]" << std::endl;

  // /* Calculate the inference error. */
  // data_.calculateStateError();
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

      utias::mrclam::utils::normaliseAngle(normalised_angle);

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

} // namespace CL::filter
