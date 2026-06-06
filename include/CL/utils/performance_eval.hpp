/**
 * @file performance_eval.hpp
 */
#pragma once

#include "CL/agent/robot.hpp"
#include "UtiasMrclam/agents/Robot.hpp"

#include <UtiasMrclam/DataHandler.hpp>
#include <memory>

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

  static void
  populateSyncedStates(const std::vector<std::unique_ptr<Robot>> &robots,
                       utias::mrclam::Handler &data);

private:
  /**
   * Find the robot whose barcode matches the one provided.
   */
  static const Robot *
  getAssociatedRobot(utias::mrclam::Robot::Barcode barcode,
                     const std::vector<std::unique_ptr<Robot>> &robots);

  static void populateSyncedStates(const Robot &robots,
                                   utias::mrclam::Robot &data);
};
} // namespace CL::utils
