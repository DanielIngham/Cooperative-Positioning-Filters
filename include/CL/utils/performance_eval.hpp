/**
 * @file performance_eval.hpp
 */
#pragma once

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

private:
};
} // namespace CL::utils
