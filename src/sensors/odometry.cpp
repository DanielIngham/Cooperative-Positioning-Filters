/**
 * @file odometry.cpp
 */

#include "CL/sensors/odometry.hpp"
#include "CL/common/types.hpp"

#include <UtiasMrclam/agents/Robot.hpp>

namespace CL::sensors {

Odometry::Odometry(
    const std::vector<utias::mrclam::Robot::Odometry> &odometry_list,
    double var_fvel, double var_avel) {

  cov_(FORWARD_VELOCITY, FORWARD_VELOCITY) = var_fvel;
  cov_(ANGULAR_VELOCITY, ANGULAR_VELOCITY) = var_avel;

  for (const auto &odometry : odometry_list) {
    data_.emplace_back(odometry.time, odometry.forward_velocity,
                       odometry.angular_velocity, cov_);
  }
}

double Odometry::timeAt(size_t index) const {
  OdomData const &odom_data{data_.at(index)};

  return odom_data.time();
}

OdomData const &Odometry::odomAt(size_t index) const { return data_.at(index); }

} // namespace CL::sensors
