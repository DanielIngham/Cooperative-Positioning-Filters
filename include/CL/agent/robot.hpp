/**
 * @file Agent.hpp
 */
#pragma once

#include "CL/agent/agent.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/ekf.hpp"
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

  template <typename FilterType = CL::filter::EKF, typename RobotType = Robot>
  static std::unique_ptr<Robot> create(const Data::Robot &data) {
    static_assert(std::is_base_of_v<filter::Filter, FilterType>,
                  "FilterType must be derived from base class filter::Filter.");

    static_assert(std::is_base_of_v<Robot, RobotType>,
                  "Type T must be derived from base class filter::Filter.");

    std::unique_ptr<Robot> robot{new RobotType{data}};

    robot->filter_ = std::make_unique<FilterType>(robot->getPrior());

    return robot;
  }

  /**
   * Uses odometry to compute the predictive density for its state, create a
   * "broadcast" its estimate over the vanet.
   * @param index The current index in the UTIAS MRCLAM synced set.
   */
  const EstimationParameters &broadcastEstimate(size_t index) override;

  /**
   * Getter for the messages (Estimates) broadcasted by other agents on the
   * VANET.
   * @param index Time index.
   * @param vanet_msgs Map containing the agent barcode (key) and estimate
   * (value).
   */
  void recieveVanetMessages(
      size_t index, std::map<unsigned short, EstimationParameters> &vanet_msgs);

  /**
   * @returns the list of the estimates produced by the agent for its entire
   * trajectory.
   */
  const std::vector<EstimationParameters> &getEstimates() const;

  const EstimationParameters &getPrior() { return estimates_.front(); }

protected:
  Robot(const Data::Robot &data);

private:
  std::unique_ptr<filter::Filter> filter_;
  std::vector<EstimationParameters> estimates_;
  const std::vector<Data::Robot::Odometry> &odometry_;
  const std::vector<Data::Robot::Measurement> &measurements_;

  double sample_period_;
};
} // namespace CL
