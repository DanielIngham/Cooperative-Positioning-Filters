/**
 * @file Inference.hpp
 */
#pragma once

#include "CL/Landmark.hpp"
#include "CL/Robot.hpp"
#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <vector>

namespace CL {

template <typename T> class Inference {
  static_assert(std::is_base_of_v<filter::Filter, T>,
                "Filter must derive from Filter abstract base class");

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

    Data::Robot::List &fleet_data{data.getRobots()};

    for (const Data::Robot &robot_data : fleet_data) {
      robots_.emplace_back(robot_data);
    }

    const Data::Landmark::List &landmarks{data.getLandmarks()};

    for (const Data::Landmark &landmark_data : landmarks) {
      landmarks_.emplace_back(landmark_data);
    }
  }

private:
  std::vector<Robot<T>> robots_;
  std::vector<Landmark> landmarks_;
};
} // namespace CL
