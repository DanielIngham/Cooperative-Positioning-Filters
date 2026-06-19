#include "CL/agent/robot.hpp"
#include "CL/agent/agent.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include "CL/sensors/meas_set.hpp"
#include "CL/sensors/odometry.hpp"
#include "CL/utils/utils.hpp"

#include <UtiasMrclam/utils/Utils.hpp>
#include <cassert>
#include <memory>

namespace CL {
Robot::Robot(const utias::mrclam::Robot &data) : Agent(data.barcode()) {

  EstimationParameters &prior{
      estimates_.emplace_back(data.id(), data.barcode())};

  /* Initial state: 3x1 Matrix. */
  prior.state_estimate(X) = data.groundtruth.states.front().x;
  prior.state_estimate(Y) = data.groundtruth.states.front().y;
  prior.state_estimate(ORIENTATION) =
      data.groundtruth.states.front().orientation;

  /* Populate the Initial information vector */
  prior.precision_matrix = prior.error_covariance.inverse();
  prior.information_vector = prior.precision_matrix * prior.state_estimate;

  odometry_ = std::make_unique<sensors::Odometry>(
      data.synced.odometry, data.forward_velocity_error.variance,
      data.angular_velocity_error.variance);

  new_measurements_ = std::make_unique<sensors::Measurements>(
      data.synced.measurements, data.range_error.variance,
      data.bearing_error.variance);
}

const EstimationParameters &Robot::broadcastEstimate(size_t index) {
  /* Assumes that the system has a known prior. */
  assert(!estimates_.empty());

  /* Check if the robot has an estimate for the current timestamp. If not
   * predict up to the timestamp.
   * NOTE: The assertion ensures that the estimates are not empty and therefore
   * this subtraction should be a safe operation.
   */
  size_t const &last_idx{estimates_.size() - 1};
  for (size_t i{last_idx}; i < index; i++) {

    double const &next_time{odometry_->timeAt(i + 1)};

    sensors::OdomData const &odometry{odometry_->odomAt(i)};

    double const dt{next_time - odometry.time()};

    assert(dt > 0);

    estimates_.emplace_back(
        filter_->prediction(odometry, estimates_.back(), dt));
  }

  return estimates_.back();
};

void Robot::recieveVanetMessages(
    size_t index, std::map<unsigned short, EstimationParameters> &vanet_msgs) {

  /* Loop through the measurements taken and perform the measurement
   * update for each robot.
   * NOTE: This operation uses the assumption that the measurements for the
   * indpendent robots/landmarks are independent of one another. */
  double const &time{odometry_->timeAt(index)};

  sensors::MeasSet const *new_curr_meas{new_measurements_->find(time)};

  if (new_curr_meas == nullptr)
    return;

  EstimationParameters &parameters{estimates_.at(index)};

  /* Loop through the list of measurements. It is assumed that the data
   * associations are known. */
  for (const sensors::MeasData &meas : *new_curr_meas) {

    auto measured_agent{vanet_msgs.find(meas.barcode())};
    if (measured_agent == vanet_msgs.end())
      continue;

    filter_->correction(parameters, measured_agent->second, meas);

    utils::normaliseAngle(parameters.state_estimate(ORIENTATION));
  }
}

const std::vector<EstimationParameters> &Robot::getEstimates() const {
  return estimates_;
}
} // namespace CL
