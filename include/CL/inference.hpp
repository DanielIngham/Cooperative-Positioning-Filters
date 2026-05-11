/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"
#include "CL/landmark.hpp"
#include "CL/robot.hpp"
#include "CL/utils/utils.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/utils/Utils.hpp>
#include <chrono>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace CL {
template <typename T> class Inference {
  static_assert(
      std::is_base_of_v<filter::Filter, T>,
      "Type passed to Inference class must be child of filter class.");

public:
  Inference() = delete;
  Inference(Inference &&) = default;
  Inference(const Inference &) = default;
  Inference &operator=(Inference &&) = default;
  Inference &operator=(const Inference &) = default;
  ~Inference() = default;

  /**
   * Populates robots and landmarks with data from the data handler.
   */
  Inference(Data::Handler &data) {
    total_datapoints_ = data.getNumberOfSyncedDatapoints();

    Data::Robot::List &fleet_data{data.getRobots()};

    for (const Data::Robot &robot_data : fleet_data) {
      robots_.emplace_back(std::make_unique<T>(data), robot_data);
    }

    const Data::Landmark::List &landmarks{data.getLandmarks()};

    for (const Data::Landmark &landmark_data : landmarks) {
      landmarks_.emplace_back(landmark_data);
    }
  }

  void compute() {
    /* Start the timer for measuring the execution time of a child filter. */
    const auto timer_start{std::chrono::high_resolution_clock::now()};

    for (size_t k{1}; k < total_datapoints_; k++) {
      std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints_
                << " %" << std::flush;

      /* Set of robot state estimates broadcasted over the VANET.
       * TODO: maybe make this an attribute of the class and turn it into a map.
       */
      std::map<unsigned short, EstimationParameters> vanet_broadcasts{};

      for (Landmark &landmark : landmarks_) {
        vanet_broadcasts[landmark.getBarcode().val()] =
            landmark.broadcastEstimate(k);
      }
      for (Robot &robot : robots_) {
        vanet_broadcasts[robot.getBarcode().val()] = robot.broadcastEstimate(k);

        // ParameterList &parameter_list{robot_parameters.at(robot.id())};
        // /* Extend the time-series by duplicating the last element in
        // place.
        //  */
        // parameter_list.push_back(parameter_list.back());
        //
        // Data::Robot::Odometry &odometry{robot.synced.odometry[k]};
        // EstimationParameters &parameters{parameter_list.back()};
        //
        // prediction(odometry, parameters);
        //
        // double normalised_angle{parameters.state_estimate(ORIENTATION)};
        // Data::Robot::normaliseAngle(normalised_angle);
        /* Update the robot state data structure. */
        // synced_states[index] = Data::Robot::State{
        //     .time = robot.groundtruth.states[k].time,
        //     .x = parameters.state_estimate(X),
        //     .y = parameters.state_estimate(Y),
        //     .orientation = normalised_angle,
        // };
      }

      for (auto &robot : robots_) {
        robot.recieveVanetMessages(k);

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
        EstimationParameters &parameters{
            robot_parameters.at(robot.id()).back()};

        for (unsigned short j{}; j < current_measurement->subjects.size();
             j++) {

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

          utils::normaliseAngle(normalised_angle);

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
  }

private:
  std::vector<Robot> robots_;
  std::vector<Landmark> landmarks_;
  size_t total_datapoints_{};
};
} // namespace CL
