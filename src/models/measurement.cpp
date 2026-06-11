#include "CL/models/measurement.hpp"
#include "CL/common/types.hpp"
#include "CL/utils/utils.hpp"
#include <algorithm>
#include <cmath>

namespace CL::Models {

Measurement::Measurement(const state_t &ego_state, const state_t &agent_state) {
  predicted_measurement_ = measurementModel(ego_state, agent_state);
  ego_jacobian_ = egoMeasurementJacobian(ego_state, agent_state);
  agent_jacobian_ = agentMeasurementJacobian(ego_state, agent_state);
  augmented_jacobian = calculateMeasurementJacobian(ego_state, agent_state);
}

measurement_t Measurement::measurementModel(const state_t &ego_state,
                                            const state_t &agent_state) {

  /* Calculate the terms */
  const double x_difference{agent_state(X) - ego_state(X)};
  const double y_difference{agent_state(Y) - ego_state(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  /* Prevent division by zero and floating point precision errors. */
  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  /* Calculate the predicted measurement based on the estimated states. */
  measurement_t predicted_measurement{
      (measurement_t() << std::sqrt((x_difference * x_difference) +
                                    (y_difference * y_difference)),
       std::atan2(y_difference, x_difference) - ego_state(ORIENTATION))
          .finished()};

  utils::normaliseAngle(predicted_measurement(BEARING));

  return predicted_measurement;
}

double Measurement::rangeMeasurementModel(const state_t &agent,
                                          const state_t &ego) {
  const double x_difference{agent(X) - ego(X)};
  const double y_difference{agent(Y) - ego(Y)};

  double range{
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2))};

  return range;
}

augmentedMeasurementJacobian_t
Measurement::calculateMeasurementJacobian(const state_t &ego_robot,
                                          const state_t &obs_agent) {

  augmentedMeasurementJacobian_t measurement_jacobian{};

  const double x_difference{obs_agent(X) - ego_robot(X)};
  const double y_difference{obs_agent(Y) - ego_robot(Y)};

  double denominator{std::hypot(x_difference, y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;

  return measurement_jacobian;
}

measurementJacobian_t
Measurement::egoMeasurementJacobian(const state_t &ego, const state_t &agent) {

  const double x_diff{agent(X) - ego(X)}, y_diff{agent(Y) - ego(Y)};
  double r{std::hypot(x_diff, y_diff)};

  static constexpr double min_r{1e-6};
  double r2{std::max(std::pow(r, 2), std::pow(min_r, 2))};

  measurementJacobian_t jacobian{};
  jacobian(RANGE, X) = -x_diff / r;
  jacobian(RANGE, Y) = -y_diff / r;
  jacobian(RANGE, ORIENTATION) = 0;

  jacobian(BEARING, X) = y_diff / r2;
  jacobian(BEARING, Y) = -x_diff / r2;
  jacobian(BEARING, ORIENTATION) = -1;

  return jacobian;
}

measurementJacobian_t
Measurement::agentMeasurementJacobian(const state_t &ego,
                                      const state_t &agent) {
  const double x_diff{agent(X) - ego(X)}, y_diff{agent(Y) - ego(Y)};
  double r{std::hypot(x_diff, y_diff)};

  static constexpr double min_r{1e-6};
  double r2{std::max(std::pow(r, 2), std::pow(min_r, 2))};

  measurementJacobian_t jacobian{};
  jacobian(RANGE, X) = x_diff / r;
  jacobian(RANGE, Y) = y_diff / r;
  jacobian(RANGE, ORIENTATION) = 0;

  jacobian(BEARING, X) = -y_diff / r2;
  jacobian(BEARING, Y) = x_diff / r2;
  jacobian(BEARING, ORIENTATION) = 0;

  return jacobian;
}

vector3D_t
Measurement::egoRangeMeasurementJacobian(const EstimationParameters &ego,
                                         const EstimationParameters &agent) {
  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double range{
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2))};

  static constexpr double min_range{1e-6};
  range = std::max(range, min_range);

  vector3D_t jacobian{};
  jacobian(X) = -x_difference / range;
  jacobian(Y) = -y_difference / range;
  jacobian(ORIENTATION) = .0;

  return jacobian;
}

vector3D_t
Measurement::agentRangeMeasurementJacobian(const EstimationParameters &ego,
                                           const EstimationParameters &agent) {
  const double x_difference{ego.state_estimate(X) - agent.state_estimate(X)};
  const double y_difference{ego.state_estimate(Y) - agent.state_estimate(Y)};
  double range =
      std::sqrt(std::pow(x_difference, 2) + std::pow(y_difference, 2));

  static constexpr double min_range{1e-6};
  range = std::max(range, min_range);

  vector3D_t jacobian;
  jacobian(X) = x_difference / range;
  jacobian(Y) = y_difference / range;
  jacobian(ORIENTATION) = .0;

  return jacobian;
}

measurement_t const &Measurement::predictedMeasurement() {
  return predicted_measurement_;
}

measurementJacobian_t const &Measurement::egoJacobian() {
  return ego_jacobian_;
}

measurementJacobian_t const &Measurement::agentJacobian() {
  return agent_jacobian_;
}

augmentedMeasurementJacobian_t const &Measurement::augmentedJacobian() {
  return augmented_jacobian;
}

} // namespace CL::Models
