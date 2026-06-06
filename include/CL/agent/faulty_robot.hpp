/**
 * @file faulty_robot.hpp
 * @author Daniel Ingham
 * @date 2026-06-06
 */
#pragma once

#include "CL/agent/robot.hpp"

namespace CL {

/**
 * A mobile robot agent in the VANET that attempts to localise its self under
 * malformed intial conditions. The agents starts off with an incorrect
 * assumption about its position.
 */
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
  /**
   * Creates a robot initialised with a malformed prior. This agent represents
   * case where an agents estimator has become inconsistent (due either
   * unaccounted cross-correlation between agent observations in the VANET,
   * badly categorised error distributions, or linearisation errors)
   * @param data Agent robot data provided by the dataset/simulation.
   */
  FaultyRobot(const utias::mrclam::Robot &data);
};
} // namespace CL
