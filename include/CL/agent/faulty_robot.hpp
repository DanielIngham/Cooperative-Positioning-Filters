/**
 * @file faulty_robot.hpp
 */
#pragma once

#include "CL/agent/robot.hpp"

namespace CL {

class FaultyRobot : public Robot {
  friend class Robot;

public:
  FaultyRobot() = delete;
  FaultyRobot(FaultyRobot &&) = default;
  FaultyRobot(const FaultyRobot &) = delete;
  FaultyRobot &operator=(FaultyRobot &&) = delete;
  FaultyRobot &operator=(const FaultyRobot &) = delete;
  ~FaultyRobot() = default;

private:
  FaultyRobot(const Data::Robot &data);
};
} // namespace CL
