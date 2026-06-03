#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

namespace CL::Models {
class Measurement {
public:
  Measurement() = delete;
  Measurement(Measurement &&) = delete;
  Measurement(const Measurement &) = delete;
  Measurement &operator=(Measurement &&) = delete;
  Measurement &operator=(const Measurement &) = delete;
  ~Measurement() = default;

  Measurement(const state_t &ego_state, const state_t &agent_state);

  /**
   * @brief Uses the non-linear measurement model to predict what the
   * measurement from the system would be given the state estimates of the ego
   * robot and measured agent.
   *
   * @param[in,out] ego_robot The estimation parameters of the ego robot.
   * @param[in] other_agent The estimation parameters of the agent that was
   * measured by the ego robot.
   *
   * @details The measusurement model for the measurement taken from ego vehicle
   * \f$i\f$ to agent \f$j\f$ used for the correction step takes the form
   * \f[ \begin{bmatrix} r_{ij}^{(t)} \\ \phi_{ij}^{(t)}\end{bmatrix} =
   * \begin{bmatrix}\sqrt{(x_j^{(t)} - x_i^{(t)})^2 + (y_j^{(t)} - y_i^{(t)})^2}
   * + q_r \\ \text{atan2}\left(\frac{y_j^{(t)}-y_i^{(t)}}{x_j^{(t)}-x_i^{(t)}
   * }\right) - \theta_i^{(t)} + q_\phi\end{bmatrix}, \f] where \f$x\f$ and
   * \f$y\f$ denote the robots coordinates; \f$\theta\f$ denotes the ego robots
   * orientation (heading); and \f$q_r\f$ and \f$q_\omega\f$ denote the Gaussian
   * distributed measurement noise.
   */
  [[nodiscard]] static measurement_t
  measurementModel(const state_t &ego_state, const state_t &agent_state);

  [[nodiscard]] static double rangeMeasurementModel(const state_t &agent,
                                                    const state_t &ego);

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
  static void calculateMeasurementJacobian(EstimationParameters &,
                                           const EstimationParameters &);

  /**
   * Calculates the Jacobian matrix of the measurement model with respect to the
   * states of the ego agent (i.e the agent that produced the measurement).
   * @param ego State of the ego vehicle.
   * @param ego State of the agent vehicle.
   */
  [[nodiscard]] static measurementJacobian_t
  egoMeasurementJacobian(const state_t &ego, const state_t &agent);

  /**
   * Calculates the Jacobian matrix of the measurement model with respect to the
   * states of the measured agent.
   * @param ego State of the ego vehicle.
   * @param ego State of the agent vehicle.
   */
  [[nodiscard]] static measurementJacobian_t
  agentMeasurementJacobian(const state_t &ego, const state_t &agent);

  [[nodiscard]] static vector3D_t
  egoRangeMeasurementJacobian(const EstimationParameters &,
                              const EstimationParameters &);

  [[nodiscard]] static vector3D_t
  agentRangeMeasurementJacobian(const EstimationParameters &,
                                const EstimationParameters &);

  /**
   * @returns the predicted measurement given the state of the ego agent and the
   * observed agent.
   */
  measurement_t predictedMeasurement();
  /**
   * @returns the Jacobian matrix of the measurement model with ego agent's
   * states.
   */
  measurementJacobian_t egoJacobian();
  /**
   * @returns the Jacobian matrix of the measurement model with observed agent's
   * states.
   */
  measurementJacobian_t agentJacobian();

private:
  measurement_t predicted_measurement_{};
  measurementJacobian_t ego_jacobian_{};
  measurementJacobian_t agent_jacobian_{};
};
} // namespace CL::Models
