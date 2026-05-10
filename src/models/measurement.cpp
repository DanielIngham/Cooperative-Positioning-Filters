#include "CL/models/measurement.hpp"
#include "CL/utils/utils.hpp"

namespace CL::Models {

/**
 * @brief Uses the non-linear measurement model to predict what the measurement
 * from the system would be given the state estimates of the ego robot and
 * measured agent.
 *
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 *
 * @details The measusurement model for the measurement taken from ego vehicle
 * \f$i\f$ to agent \f$j\f$ used for the correction step takes the form
 * \f[ \begin{bmatrix} r_{ij}^{(t)} \\ \phi_{ij}^{(t)}\end{bmatrix} =
 * \begin{bmatrix}\sqrt{(x_j^{(t)} - x_i^{(t)})^2 + (y_j^{(t)} - y_i^{(t)})^2} +
 * q_r \\ \text{atan2}\left(\frac{y_j^{(t)}-y_i^{(t)}}{x_j^{(t)}-x_i^{(t)}
 * }\right) - \theta_i^{(t)} + q_\phi\end{bmatrix}, \f] where \f$x\f$ and
 * \f$y\f$ denote the robots coordinates; \f$\theta\f$ denotes the ego robots
 * orientation (heading); and \f$q_r\f$ and \f$q_\omega\f$ denote the Gaussian
 * distributed measurement noise (See
 * Filter::EstimationParameters.measurement_noise).
 */
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

/**
 * @brief Jacobian of the measurement model evaluated in terms of the systems
 * states: x,y, and heading.
 *
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the agent that was
 * measured by the ego robot.
 *
 * @details The formula used for the calculation of the Jacobian of the
 * measurement matrix between ego vehicle \f$i\f$ and measured agent
 * \f$j\f$ take the form
 * \f[ H = \begin{bmatrix} \frac{-\Delta x}{d} & \frac{-\Delta y}{d} & 0 &
 * \frac{\Delta x}{d} & \frac{\Delta y}{d} & 0\\ \frac{\Delta y}{d^2} &
 * \frac{-\Delta x}{d^2} & -1 & \frac{-\Delta y}{d^2} & \frac{\Delta x}{d^2}
 * & 0\end{bmatrix} \f] where \f$\Delta x = x_j - x_i\f$; \f$\Delta y = y_j
 * - y_i\f$; and \f$\Delta d = \sqrt{\Delta x^2 + \Delta y^2}\f$.
 */
void Measurement::calculateMeasurementJacobian(
    EstimationParameters &ego_robot, const EstimationParameters &other_agent) {

  const double x_difference{other_agent.state_estimate(X) -
                            ego_robot.state_estimate(X)};

  const double y_difference{other_agent.state_estimate(Y) -
                            ego_robot.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  ego_robot.measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;
}

measurementJacobian_t
Measurement::egoMeasurementJacobian(const EstimationParameters &ego,
                                    const EstimationParameters &agent) {

  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  measurementJacobian_t jacobian;
  jacobian(0, 0) = -x_difference / denominator;
  jacobian(0, 1) = -y_difference / denominator;
  jacobian(0, 2) = 0;

  jacobian(1, 0) = y_difference / (denominator * denominator);
  jacobian(1, 1) = -x_difference / (denominator * denominator);
  jacobian(1, 2) = -1;
  return jacobian;
}

measurementJacobian_t
Measurement::agentMeasurementJacobian(const EstimationParameters &ego,
                                      const EstimationParameters &agent) {
  const double x_difference{agent.state_estimate(X) - ego.state_estimate(X)};

  const double y_difference{agent.state_estimate(Y) - ego.state_estimate(Y)};

  double denominator{
      std::sqrt(x_difference * x_difference + y_difference * y_difference)};

  static constexpr double min_distance{1e-6};
  if (denominator < min_distance)
    denominator = min_distance;

  measurementJacobian_t jacobian;
  jacobian(0, 0) = x_difference / denominator;
  jacobian(0, 1) = y_difference / denominator;
  jacobian(0, 2) = 0;

  jacobian(1, 0) = -y_difference / (denominator * denominator);
  jacobian(1, 1) = x_difference / (denominator * denominator);
  jacobian(1, 2) = 0;

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

  vector3D_t jacobian;
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

} // namespace CL::Models
