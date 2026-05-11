#include "CL/robot.hpp"
#include "CL/agent.hpp"

#include "CL/common/estimation_parameters.hpp"
#include "CL/utils/utils.hpp"
#include <UtiasMrclam/utils/Utils.hpp>

#include <cassert>
#include <iostream>

namespace CL {

Robot::Robot(std::unique_ptr<filter::Filter> filter_ptr,
             const Data::Robot &data)
    : Agent(data.barcode()), filter_{std::move(filter_ptr)},
      odometry_{data.synced.odometry}, measurements_{data.synced.measurements} {

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
  /* Assumes that the system has a known prior. */
  assert(!estimates_.empty());

  /* Check if the robot has an estimate for the current timestamp. If not
   * predict up to the timestamp.
   * NOTE: The assertion ensures that the estimates are not empty and therefore
   * this subtraction should be a safe operation.
   */
  for (size_t i{estimates_.size() - 1}; i <= index; i++) {

    /* Extend the time-series by duplicating the last element in place. */
    auto &current_estimate{estimates_.emplace_back(estimates_.back())};

    const Data::Robot::Odometry &current_odometry{odometry_.at(i)};

    filter_->prediction(current_odometry, current_estimate);
  }

  return estimates_.back();
};

void Robot::recieveVanetMessages(
    size_t index, std::map<unsigned short, EstimationParameters> &vanet_msgs) {

  /* Loop through the measurements taken and perform the measurement
   * update for each robot.
   * NOTE: This operation uses the assumption that the measurements fo the
   * indpendent robots/landmarks are independent of one another.
   */
  const double time{odometry_.at(index).time};
  const Data::Robot::Measurement *current_measurement{
      utias::mrclam::utils::getMeasurement(measurements_, time)};

  if (current_measurement == nullptr) {
    return;
  }

  EstimationParameters &parameters{estimates_.at(index)};

  std::cerr << "Before\n" << parameters.state_estimate << std::endl;
  /* Loop through the list of measurements. It is assumed that the data
   * associations are known. */
  for (unsigned short i{}; i < current_measurement->subjects.size(); i++) {

    /* Find the subject for whom the barcode belongs to and check if they are on
     * the VANET. */
    const Data::Agent::Barcode &barcode{current_measurement->subjects.at(i)};
    auto measured_agent{vanet_msgs.find(barcode.val())};

    if (measured_agent == vanet_msgs.end()) {
      continue;
    }

    parameters.measurement[RANGE] = current_measurement->ranges.at(i);
    parameters.measurement[BEARING] = current_measurement->bearings.at(i);

    filter_->correction(parameters, measured_agent->second);

    double &normalised_angle{parameters.state_estimate(ORIENTATION)};
    utils::normaliseAngle(normalised_angle);
  }
  std::cerr << "After\n" << parameters.state_estimate << std::endl;
}

const std::vector<EstimationParameters> &Robot::getEstimates() const {
  return estimates_;
}
} // namespace CL
