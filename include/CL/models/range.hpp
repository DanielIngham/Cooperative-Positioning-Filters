#include "CL/models/measurement.hpp"

namespace CL::Models {

/**
 * Implements the functionality of the range only measurement model for
 * filtering. This includes:
 * - Computing the predicted measurement given the state of the ego vehicle and
 * the observed agent.
 * - Computing the corresponding Jacobian matrix.
 */
class Range : public Measurement {
public:
  Range() = delete;
  Range(Range &&) = default;
  Range(const Range &) = default;
  Range &operator=(Range &&) = default;
  Range &operator=(const Range &) = default;
  ~Range() = default;

  Range(const EstimationParameters &ego, const EstimationParameters agent);

private:
  [[nodiscard]] static Eigen::Matrix<double, 1, 3>
  egoRangeMeasurementJacobian(const EstimationParameters &ego,
                              const EstimationParameters &agent);

  [[nodiscard]] double model(const state_t &, const state_t &);
};
} // namespace CL::Models
