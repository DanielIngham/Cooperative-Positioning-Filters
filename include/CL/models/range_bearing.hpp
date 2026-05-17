/**
 * @file range_bearing.hpp
 */
#pragma once
#include "CL/common/types.hpp"
#include "CL/models/measurement.hpp"

namespace CL::Models {

/**
 * Implements the functionality of the range and bearing measurement model for
 * filtering. This includes:
 * - Computing the predicted measurement given the state of the ego vehicle and
 * the observed agent.
 * - Computing the corresponding Jacobian matrix.
 */
class RangeBearing : public Measurement {
public:
  RangeBearing() = delete;
  RangeBearing(RangeBearing &&) = default;
  RangeBearing(const RangeBearing &) = default;
  RangeBearing &operator=(RangeBearing &&) = default;
  RangeBearing &operator=(const RangeBearing &) = default;
  ~RangeBearing() = default;

  RangeBearing(const EstimationParameters &ego,
               const EstimationParameters agent);

private:
  [[nodiscard]] measurement_t model(const state_t &ego,
                                    const state_t &agent) override;

  [[nodiscard]] augmentedMeasurementJacobian_t
  calculateMeasurementJacobian(const EstimationParameters &ego,
                               const EstimationParameters &agent) const;
};
} // namespace CL::Models
