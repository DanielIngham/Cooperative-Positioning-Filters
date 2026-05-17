#include "CL/models/measurement.hpp"

namespace CL::Models {
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
  measurement_t model(const state_t &ego, const state_t &agent) override;

  void calculateMeasurementJacobian(const EstimationParameters &ego,
                                    const EstimationParameters &agent);
};
} // namespace CL::Models
