/**
 * @file Agent.hpp
 */
#pragma once

#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Robot.hpp>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace CL {

class Robot {

public:
  Robot() = delete;
  Robot(Robot &&) = default;
  Robot(const Robot &) = delete;
  Robot &operator=(Robot &&) = delete;
  Robot &operator=(const Robot &) = delete;
  ~Robot() = default;

  Robot(std::unique_ptr<filter::Filter> filter_ptr, Data::Robot &data);

  const EstimationParameters &broadcastEstimate(size_t index);

private:
  std::unique_ptr<filter::Filter> filter_;
  std::vector<EstimationParameters> estimates_;
  std::vector<Data::Robot::Odometry> &odometry_;
  std::vector<Data::Robot::Measurement> &measurements_;
};
} // namespace CL
