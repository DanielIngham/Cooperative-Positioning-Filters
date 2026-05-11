/**
 * @file performance_eval.hpp
 */
#pragma once

#include "CL/robot.hpp"
#include "UtiasMrclam/agents/Robot.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <vector>

namespace CL::utils {
class PerformanceEvaluator {
public:
  PerformanceEvaluator() = default;
  PerformanceEvaluator(PerformanceEvaluator &&) = default;
  PerformanceEvaluator(const PerformanceEvaluator &) = default;
  PerformanceEvaluator &operator=(PerformanceEvaluator &&) = default;
  PerformanceEvaluator &operator=(const PerformanceEvaluator &) = default;
  ~PerformanceEvaluator() = default;

  void writeInnovation();
  void writeNormalisedInnovation();
  void writeNEES();

  static void populateSyncedStates(const std::vector<Robot> &robots,
                                   Data::Handler &data);

private:
  /**
   * Find the robot whose barcode matches the one provided.
   */
  static const Robot *getAssociatedRobot(Data::Robot::Barcode barcode,
                                         const std::vector<Robot> &robots);

  static void populateSyncedStates(const Robot &robots, Data::Robot &data);
};
} // namespace CL::utils
