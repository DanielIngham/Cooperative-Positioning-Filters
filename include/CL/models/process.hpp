#pragma once

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

namespace CL::Models {
class Process {
public:
  Process() = default;
  Process(Process &&) = delete;
  Process(const Process &) = delete;
  Process &operator=(Process &&) = delete;
  Process &operator=(const Process &) = delete;
  ~Process() = default;

  static void motionModel(const Data::Robot::Odometry &, state_t &,
                          const double);

  static void calculateMotionJacobian(const Data::Robot::Odometry &,
                                      EstimationParameters &, const double);

  static void calculateProcessJacobian(EstimationParameters &, const double);

private:
};

} // namespace CL::Models
