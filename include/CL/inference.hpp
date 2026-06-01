/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/agent/adversarial_landmark.hpp"
#include "CL/agent/faulty_robot.hpp"
#include "CL/agent/landmark.hpp"
#include "CL/agent/robot.hpp"
#include "CL/common/config.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/utils/Utils.hpp>
#include <chrono>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace CL {
template <typename FilterType> class Inference {
  static_assert(
      std::is_base_of_v<filter::Filter, FilterType>,
      "Type passed to Inference class must be child of filter class.");

public:
  Inference() = delete;
  Inference(Inference &&) = default;
  Inference(const Inference &) = default;
  Inference &operator=(Inference &&) = default;
  Inference &operator=(const Inference &) = default;
  ~Inference() = default;

  /**
   * Populates robots and landmarks with data from the data handler.
   * @param data Dataset.
   */
  Inference(Data::Handler &data, const Config &config = {}) {

    total_timesteps_ = data.getNumberOfSyncedDatapoints();

    Data::Robot::List &fleet_data{data.getRobots()};

    size_t i{};
    for (const Data::Robot &robot_data : fleet_data) {
      /* Populate the VANET first with cooperative robots, then with faulty
       * robots. */
      if (i++ < config.robots.cooperative)
        robots_.push_back(Robot::create<FilterType, Robot>(robot_data));
      else
        robots_.push_back(Robot::create<FilterType, FaultyRobot>(robot_data));
    }

    const Data::Landmark::List &landmarks{data.getLandmarks()};

    size_t j{};
    for (const Data::Landmark &landmark_data : landmarks) {
      // if (j++ % 4 == 0) {
      if (j++ < config.landmarks.adversarial) {
        landmarks_.emplace_back(
            std::make_unique<AdversarialLandmark>(landmark_data));
      } else
        landmarks_.emplace_back(std::make_unique<Landmark>(landmark_data));
    }
  }

  /**
   * @returns a vector containing the robot agents in the VANET.
   */
  const std::vector<std::unique_ptr<Robot>> &getRobots() { return robots_; }

  /**
   * Facilitates the distributed state estimates in the VANET through the
   * passing of state-estimates between agents in the VANET.
   */
  void compute() {
    /* Start the timer for measuring the execution time of a child filter. */
    const auto timer_start{std::chrono::high_resolution_clock::now()};

    for (size_t k{1}; k < total_timesteps_; ++k) {
      std::cout << "\rPerforming Inference: " << k * 100 / total_timesteps_
                << " %" << std::flush;

      receiveVanetMessages(k);
      distributeVanetMessages(k);
    }

    const auto timer_end{std::chrono::high_resolution_clock::now()};
    const auto execution_duration{
        std::chrono::duration_cast<std::chrono::milliseconds>(timer_end -
                                                              timer_start)};
    std::cout << '\r' << std::flush;
    std::cout << "\033[1;32mInference Complete:\033[0m ["
              << execution_duration.count() << " ms]" << std::endl;
  }

private:
  /**
   * List of robot agents in the VANET.
   */
  std::vector<std::unique_ptr<Robot>> robots_;

  /**
   * List of landmark agents in the VANET.
   */
  std::vector<std::unique_ptr<Landmark>> landmarks_;

  /**
   * Total number of time steps to loop through in the dataset.
   */
  size_t total_timesteps_{};

  /**
   * Map of messages broadcasted over the VANET.
   */
  std::map<unsigned short, EstimationParameters> vanet_broadcasts_{};

  /**
   * Collects all the messages broadcasted over the VANET to be passed to the
   * agents. This simulates broadcasting over a network.
   * @param index Time index in the dataset.
   */
  void receiveVanetMessages(size_t index) {
    for (auto &landmark : landmarks_) {
      vanet_broadcasts_[landmark->getBarcode().val()] =
          landmark->broadcastEstimate(index);
    }

    for (auto &robot_ptr : robots_) {
      vanet_broadcasts_[robot_ptr->getBarcode().val()] =
          robot_ptr->broadcastEstimate(index);
    }
  }

  /**
   * Passes the broadcasted messages to each of the robots in the Vanet.
   * @param index Time index in the dataset.
   * @note Since Landmarks aren't actively estimating their position, they don't
   * need to recieve the broadcasted estimates.
   */
  void distributeVanetMessages(size_t index) {
    for (auto &robot : robots_)
      robot->recieveVanetMessages(index, vanet_broadcasts_);
  }
};
} // namespace CL
