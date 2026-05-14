/**
 * @file Agent.hpp
 */
#pragma once

#include "CL/agent.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Robot.hpp>
#include <cassert>
#include <memory>
#include <type_traits>

namespace CL {

class Robot : public Agent {

public:
  Robot() = delete;
  Robot(Robot &&) = default;
  Robot(const Robot &) = delete;
  Robot &operator=(Robot &&) = delete;
  Robot &operator=(const Robot &) = delete;
  ~Robot() = default;

  template <typename T> static Robot create(const Data::Robot &data) {
    static_assert(std::is_base_of_v<filter::Filter, T>,
                  "Type T must be derived from base class filter::Filter.");
    Robot robot{data};
    robot.filter_ = std::make_unique<T>(robot.estimates_.front());

    return robot;
  }

  /**
   * Uses odometry to compute the predictive density for its state, create a
   * "broadcast" its estimate over the vanet.
   * @param index The current index in the UTIAS MRCLAM synced set.
   */
  const EstimationParameters &broadcastEstimate(size_t index) override;
  void recieveVanetMessages(
      size_t index, std::map<unsigned short, EstimationParameters> &vanet_msgs);

  const std::vector<EstimationParameters> &getEstimates() const;

private:
  Robot(const Data::Robot &data);

  std::unique_ptr<filter::Filter> filter_;
  std::vector<EstimationParameters> estimates_;
  const std::vector<Data::Robot::Odometry> &odometry_;
  const std::vector<Data::Robot::Measurement> &measurements_;

  double sample_period_;
};
} // namespace CL
