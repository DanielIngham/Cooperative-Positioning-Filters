#include "CL/models/range.hpp"
#include "CL/common/types.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace CL::Models {
Range::Range(const EstimationParameters &ego,
             const EstimationParameters agent) {
  ego_jacobian_ = egoRangeMeasurementJacobian(ego, agent);

  agent_jacobian_ = -ego_jacobian_;

  measurement_jacobian_.resize(1, 6);
  measurement_jacobian_.block<1, 3>(0, 0) = ego_jacobian_;
  measurement_jacobian_.block<1, 3>(0, 3) = agent_jacobian_;

  predicted_measurement_ = Eigen::MatrixXd::Constant(
      1, 1, model(ego.state_estimate, agent.state_estimate));
}

Eigen::Matrix<double, 1, 3>
Range::egoRangeMeasurementJacobian(const EstimationParameters &ego,
                                   const EstimationParameters &agent) {

  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};
  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double range{std::hypot(x_difference, y_difference)};

  static constexpr double min_range{1e-6};
  range = std::max(range, min_range);

  Eigen::Matrix<double, 1, 3> jacobian{};
  jacobian(0, X) = -x_difference / range;
  jacobian(0, Y) = -y_difference / range;
  jacobian(0, ORIENTATION) = .0;

  return jacobian;
}

double Range::model(const state_t &agent, const state_t &ego) {
  const double x_difference{agent(X) - ego(X)};
  const double y_difference{agent(Y) - ego(Y)};

  double range{std::hypot(x_difference, y_difference)};

  return range;
}
} // namespace CL::Models
