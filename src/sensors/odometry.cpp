/**
 * @file odometry.cpp
 */

#include "CL/sensors/odometry.hpp"
#include "CL/common/types.hpp"

#include <UtiasMrclam/agents/Robot.hpp>
#include <vector>

namespace CL::sensors {
Odometry::Odometry(const std::vector<Data::Robot::Odometry> &odometry_list,
                   double forw_vel_var, double ang_vel_var) {
  for (const auto &odometry : odometry_list) {
  }
  time_ = odometry.time;

  input_(FORWARD_VELOCITY) = odometry_list.forward_velocity;
  input_(ANGULAR_VELOCITY) = odometry_list.angular_velocity;

  cov_(FORWARD_VELOCITY, FORWARD_VELOCITY) = forw_vel_var;
  cov_(ANGULAR_VELOCITY, ANGULAR_VELOCITY) = ang_vel_var;
}

} // namespace CL::sensors
