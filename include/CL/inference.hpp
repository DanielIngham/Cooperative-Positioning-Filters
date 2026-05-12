/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"
#include "CL/landmark.hpp"
#include "CL/robot.hpp"

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
      robots_.emplace_back(std::make_unique<T>(data), robot_data);
    }

    const Data::Landmark::List &landmarks{data.getLandmarks()};

    for (const Data::Landmark &landmark_data : landmarks) {
      landmarks_.emplace_back(landmark_data);
    }
  }

  const std::vector<Robot> &getRobots() { return robots_; }

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
  std::vector<Robot> robots_;
  std::vector<Landmark> landmarks_;
  size_t total_datapoints_{};

  std::map<unsigned short, EstimationParameters> vanet_broadcasts_{};

  void receiveVanetMessages(size_t index) {
    for (Landmark &landmark : landmarks_) {
      vanet_broadcasts_[landmark.getBarcode().val()] =
          landmark.broadcastEstimate(index);
    }

    for (Robot &robot : robots_) {
      vanet_broadcasts_[robot.getBarcode().val()] =
          robot.broadcastEstimate(index);
    }
  }

  void distributeVanetMessages(size_t index) {
    for (Robot &robot : robots_)
      robot.recieveVanetMessages(index, vanet_broadcasts_);
  }
};
} // namespace CL
