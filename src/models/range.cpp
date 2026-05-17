#include "CL/models/range.hpp"

namespace CL::Models {
Range::Range(const EstimationParameters &ego,
             const EstimationParameters agent) {
  measurement_jacobian_.leftCols<3>() = egoRangeMeasurementJacobian(ego, agent);
  measurement_jacobian_.rightCols<3>() =
      agentRangeMeasurementJacobian(ego, agent);
  predicted_measurement_ =
      rangeMeasurementModel(ego.state_estimate, agent.state_estimate);
}

vector3D_t
Range::egoRangeMeasurementJacobian(const EstimationParameters &ego,
                                   const EstimationParameters &agent) {
  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double range{
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2))};

  static constexpr double min_range{1e-6};
  range = std::max(range, min_range);

  vector3D_t jacobian;
  jacobian(X) = -x_difference / range;
  jacobian(Y) = -y_difference / range;
  jacobian(ORIENTATION) = .0;

  return jacobian;
}

vector3D_t
Range::agentRangeMeasurementJacobian(const EstimationParameters &ego,
                                     const EstimationParameters &agent) {
  const double x_difference{ego.state_estimate(X) - agent.state_estimate(X)};
  const double y_difference{ego.state_estimate(Y) - agent.state_estimate(Y)};
  double range{
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2))};

  static constexpr double min_range{1e-6};
  range = std::max(range, min_range);

  vector3D_t jacobian;
  jacobian(X) = x_difference / range;
  jacobian(Y) = y_difference / range;
  jacobian(ORIENTATION) = .0;

  return jacobian;
}

measurement_t Range::rangeMeasurementModel(const state_t &agent,
                                           const state_t &ego) {
  const double x_difference{agent(X) - ego(X)};
  const double y_difference{agent(Y) - ego(Y)};

  double range{
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2))};

  return range;
}
} // namespace CL::Models
