#include "CL/models/measurement.hpp"

namespace CL::Models {
class RangeBearing : public Measurement {
public:
  RangeBearing() = default;
  RangeBearing(RangeBearing &&) = default;
  RangeBearing(const RangeBearing &) = default;
  RangeBearing &operator=(RangeBearing &&) = default;
  RangeBearing &operator=(const RangeBearing &) = default;
  ~RangeBearing() = default;

  static measurement_t measurementModel(const state_t &, const state_t &);

private:
};
} // namespace CL::Models
