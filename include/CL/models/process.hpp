#pragma once

#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>

#include "CL/common/estimation_parameters.hpp"
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

  Process(const Data::Robot::Odometry &odometry,
          EstimationParameters &parameters, double sample_period);

  [[nodiscard]] motionJacobian_t getMotionJacobian() const;
  [[nodiscard]] processJacobian_t getProcessJacobian() const;
  [[nodiscard]] state_t getPredictedState() const;

private:
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
   * @brief Estimated robot state - x coordinate [m], y-coordinate [m], and
   * orientation (heading) [rad]: 3x1 matrix.
   * @details The state vector of the robot take the form \f[\begin{bmatrix} x
   * & y & \theta \end{bmatrix}^\top, \f] where \f$x\f$ and \f$y\f$ denotes
   * the robots 2D coordinates; and \f$\theta \f$ denotes the robots heading.
   */
  state_t state_estimate_{state_t::Zero()};

  [[nodiscard]] state_t motionModel(const Data::Robot::Odometry &, state_t &,
                                    const double);

  motionJacobian_t calculateMotionJacobian(const Data::Robot::Odometry &,
                                           EstimationParameters &,
                                           const double);

  processJacobian_t calculateProcessJacobian(const EstimationParameters &,
                                             const double);
};

} // namespace CL::Models
