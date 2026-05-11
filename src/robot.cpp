#include "CL/robot.hpp"
#include "CL/agent.hpp"

namespace CL {

Robot::Robot(std::unique_ptr<filter::Filter> filter_ptr, Data::Robot &data)
    : Agent(data.barcode()), filter_{std::move(filter_ptr)},
      odometry_{data.synced.odometry}, measurements_{data.synced.measurements},
      synced_states_{data.synced.states} {
  synced_states_.front() = data.groundtruth.states.front();

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
  parameters.measurement_noise(BEARING, BEARING) = data.bearing_error.variance;

  /* Populate the Initial information vector */
  parameters.information_vector =
      parameters.precision_matrix * parameters.state_estimate;
}

const EstimationParameters &Robot::broadcastEstimate(size_t index) {
  // TODO: Check if the robot has an estimate for the current timestamp. If not
  // predict up to the timestamp.
  // assert(index == estimates_.size() - 1);

  /* Extend the time-series by duplicating the last element in place. */
  auto &current_estimate{estimates_.emplace_back(estimates_.back())};

  const Data::Robot::Odometry &current_odometry{odometry_.at(index)};

  filter_->prediction(current_odometry, current_estimate);

  return current_estimate;
};

const EstimationParameters &Robot::recieveVanetMessages(size_t index) {}

void Robot::updateSyncedStates(size_t index,
                               const EstimationParameters &estimate) {}
} // namespace CL
