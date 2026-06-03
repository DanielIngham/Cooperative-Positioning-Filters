#pragma once

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>

#include "CL/common/types.hpp"

namespace CL::Models {
class Process {
public:
  Process() = delete;
  Process(Process &&) = delete;
  Process(const Process &) = delete;
  Process &operator=(Process &&) = delete;
  Process &operator=(const Process &) = delete;
  ~Process() = default;

  /**
   * Process model constructor.
   * @param[in] odometry The prior inputs into the system comprising a forward
   * and angular velocity.
   * @param[in] state The state estimate at the previous time index.
   * @param[in] sample_period The period between odometry measurements.
   */
  Process(const Data::Robot::Odometry &odometry, const state_t &state,
          const double sample_period);

  /**
   * @returns Mean of the predictive density.
   */
  state_t predictedState();
  /**
   * @returns Jacobian of the process model with respect to the system states.
   */
  motionJacobian_t motionJacobian();
  /**
   * @returns Jacobian of the process model with respect to the system inputs.
   */
  processJacobian_t processJacobian();

private:
  /**
   * Mean of the predictive density.
   */
  state_t predicted_state_;
  /**
   * @brief Jacobian matrix of the motion model evaluated in terms of the
   * systems states: x, y and orientation.
   */
  motionJacobian_t motion_jacobian_{motionJacobian_t::Zero()};

  /**
   * @brief Jacobian matrix of the motion model evaluated in terms of the
   * process inputs.
   */
  processJacobian_t process_jacobian_{processJacobian_t::Zero()};

  /**
   * @brief The unicycle motion model used to perform motion predictions.
   *
   * @param[in] odometry The prior inputs into the system comprising a forward
   * and angular velocity.
   * @param[in] state The state estimate at the previous time index.
   * @param[in] sample_period The period between odometry measurements.
   *
   * @details The motion model used for robot takes the form:
   * \f[\begin{bmatrix} x_i^{(t+1)} \\  y_i^{(t+1)}
   * \\ \theta_i^{(t+1)}\end{bmatrix} = \begin{bmatrix} x_i^{(t)} +
   * \tilde{v}_i^{(t)}\Delta t\cos(\theta_i^{(t)}) \\ y_i^{(t)} +
   * \tilde{v}_i^{(t)}\Delta t\sin(\theta_i^{(t)})\\ \theta_i^{(t)} +
   * \tilde{\omega}_i^{(t)}\Delta t. \end{bmatrix}, \f] where \f$i\f$ denotes
   * the robots ID; \f$t\f$ denotes the current timestep; \f$\Delta t\f$ denotes
   * the sample period; \f$x\f$ and \f$y\f$ are the robots coordinates;
   * \f$\theta\f$ denotes the robots heading (orientation); \f$\tilde{v}_i\f$
   * denotes the forward velocity input; and \f$\tilde{\omega}\f$ denotes the
   * angular velocity input. Both \f$\tilde{v}_i\f$ and \f$\tilde{\omega}_i\f$
   * are normally distributed random variables \f$\mathcal{N}(0,w)\f$ (see
   * Filter::EstimationParameters::process_noise).
   */
  state_t motionModel(const Data::Robot::Odometry &odometry,
                      const state_t &state, const double sample_period);

  /**
   * @brief Calculates the Jacobian matrix of the unicycle motion model in terms
   * of the systems states (x,y,orientation).
   *
   * @param[in] odometry The prior inputs into the system comprising a forward
   * and angular velocity.
   * @param[in] state The state estimate at the previous time index.
   * @param[in] sample_period The period between odometry measurements.
   *
   * @details The formula used for the calculation of the motion model
   * Jacobian takes the form:
   * \f[ F = \begin{bmatrix} 1 & 0 & -\tilde{v}\Delta t \sin(\theta) \\ 0 & 1
   * & \tilde{v} \Delta t \cos(\theta) \\ 0 & 0 & 1 \end{bmatrix}, \f]
   * where \f$\theta\f$ denotes the heading (orientation) of the ego vehicle;
   * and \f$\tilde{v}\f$ denotes the forward velocity. The forward velocity is
   * a random variable with Gaussian distributed noise \f$\mathcal{N}(0,w)\f$
   * , where \f$w\f$ is defined by the covariance matrix
   * Filter::EstimationParameters.measurement_noise. See Filter::motionModel for
   * information on the motion model from which this was derived.
   */
  motionJacobian_t
  calculateMotionJacobian(const Data::Robot::Odometry &odometry,
                          const state_t &state, const double sample_period);

  /**
   * @brief Calculates the Jacobian matrix of the motion model evaluated in
   * terms of the process inputs.
   *
   * @param[in] state The state estimate at the previous time index.
   * @param[in] sample_period The period between odometry measurements.
   *
   * @details The formula used for the calculation of the process noise
   * Jacobian takes the form
   * \f[L = \begin{bmatrix}\Delta t \cos(\theta) & 0 \\ \Delta t \sin(\theta)
   * & 0 \\ 0 & \Delta t \end{bmatrix}, \f] where \f$\Delta t\f$ denotes the
   * sample period; and \f$\theta\f$ denotes the heading (orientation) of the
   * ego robot. measurement_noise. See EKF::prediction for information on the
   * motion model from which this was derived.
   */
  processJacobian_t calculateProcessJacobian(const state_t &state,
                                             const double sample_period);
};

} // namespace CL::Models
