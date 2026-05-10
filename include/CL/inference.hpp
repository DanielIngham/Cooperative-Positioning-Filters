/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/Landmark.hpp"
#include "CL/Robot.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <chrono>
#include <iostream>
#include <vector>

namespace CL {

class Inference {
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
      robots_.emplace_back(robot_data);
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
      std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints
                << " %" << std::flush;

      for (auto &robot : robots_) {
        // robot.publishDataOverNetwork()?
        // create a collection of Estimation parameters based on this
        // then share that to each robot as part of the update step.
        // Kind of simulating the idea that each agent
        EstimationParameters;
        ParameterList &parameter_list{robot_parameters.at(robot.id())};
        /* Extend the time-series by duplicating the last element in place.
         */
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
  }

private:
  std::vector<Robot<T>> robots_;
  std::vector<Landmark> landmarks_;
  size_t total_datapoints_{};
};
} // namespace CL
