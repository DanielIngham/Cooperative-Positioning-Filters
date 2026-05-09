#pragma once

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>

#include "CLFilters/common/estimation_parameters.hpp"
#include "CLFilters/common/types.hpp"

namespace Filters::Models {
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

} // namespace Filters::Models
