/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/agent/adversarial_landmark.hpp"
#include "CL/agent/landmark.hpp"
#include "CL/agent/robot.hpp"
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
template <typename T> class Inference {
  static_assert(
      std::is_base_of_v<filter::Filter, T>,
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
   */
  Inference(Data::Handler &data) {

    total_datapoints_ = data.getNumberOfSyncedDatapoints();

    Data::Robot::List &fleet_data{data.getRobots()};

    for (const Data::Robot &robot_data : fleet_data) {
      robots_.push_back(Robot::create<T, Robot>(robot_data));
    }

    const Data::Landmark::List &landmarks{data.getLandmarks()};

    size_t landmark_number{};
    for (const Data::Landmark &landmark_data : landmarks) {
      if (landmark_number++ % 4 == 0) {
        landmarks_.emplace_back(
            std::make_unique<AdversarialLandmark>(landmark_data));
        std::cerr << "Add AdversarialLandmark" << std::endl;
      } else
        landmarks_.emplace_back(std::make_unique<Landmark>(landmark_data));
    }
  }

  const std::vector<std::unique_ptr<Robot>> &getRobots() { return robots_; }

  void compute() {
    /* Start the timer for measuring the execution time of a child filter. */
    const auto timer_start{std::chrono::high_resolution_clock::now()};

    for (size_t k{1}; k < total_datapoints_; ++k) {
      std::cout << "\rPerforming Inference: " << k * 100 / total_datapoints_
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
  std::vector<std::unique_ptr<Robot>> robots_;
  std::vector<std::unique_ptr<Landmark>> landmarks_;
  size_t total_datapoints_{};

  std::map<unsigned short, EstimationParameters> vanet_broadcasts_{};

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

  void distributeVanetMessages(size_t index) {
    for (auto &robot : robots_)
      robot->recieveVanetMessages(index, vanet_broadcasts_);
  }
};
} // namespace CL
