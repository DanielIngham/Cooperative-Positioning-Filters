/**
 * @file Agent.hpp
 */
#pragma once

#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Robot.hpp>

namespace CL {

template <typename T> class Robot {
  static_assert(std::is_base_of_v<filter::Filter, T>,
                "Filter must derive from Filter abstract base class");

public:
  Robot() = delete;
  Robot(Robot &&) = default;
  Robot(const Robot &) = default;
  Robot &operator=(Robot &&) = default;
  Robot &operator=(const Robot &) = default;
  ~Robot() = default;

  Robot(Data::Robot &data) {

    /* TODO: REMOVE later.
     * We want to stop saving data into the data handler down the line, but for
     * now it must remain for testing. */
    auto parameters{estimates_.emplace_back(data.id(), data.barcode())};
    /* Initial state: 3x1 Matrix. */
    parameters.state_estimate(X) = data.groundtruth.states.front().x;
    parameters.state_estimate(Y) = data.groundtruth.states.front().y;
    parameters.state_estimate(ORIENTATION) =
        data.groundtruth.states.front().orientation;

    /* Populate odometry error covariance matrix: 2x2 matrix. */
    parameters.process_noise(FORWARD_VELOCITY, FORWARD_VELOCITY) =
        data.forward_velocity_error.variance;
    parameters.process_noise(ANGULAR_VELOCITY, ANGULAR_VELOCITY) =
        data.angular_velocity_error.variance;

    /* Populate measurement error covariance matrix: 2x2 matrix. */
    parameters.measurement_noise(RANGE, RANGE) = data.range_error.variance;
    parameters.measurement_noise(BEARING, BEARING) =
        data.bearing_error.variance;

    /* Populate the Initial information vector */
    parameters.information_vector =
        parameters.precision_matrix * parameters.state_estimate;
  }

private:
  T filter_;
  std::vector<EstimationParameters> estimates_;
  std::vector<Data::Robot::Odometry> odometry_;
  std::vector<Data::Robot::Measurement> measurements_;
};
} // namespace CL
