#pragma once

#include "CL/common/types.hpp"

#include <UtiasMrclam/DataHandler.hpp>

namespace CL {

/**
 * @struct EstimationParameters
 * @brief Houses the parameters required for performing Bayesian filtering.
 */
struct EstimationParameters {
  utias::mrclam::Agent::ID id;
  utias::mrclam::Agent::Barcode barcode;

  /**
   * @brief Estimated robot state - x coordinate [m], y-coordinate [m], and
   * orientation (heading) [rad]: 3x1 matrix.
   * @details The state vector of the robot take the form \f[\begin{bmatrix} x
   * & y & \theta \end{bmatrix}^\top, \f] where \f$x\f$ and \f$y\f$ denotes
   * the robots 2D coordinates; and \f$\theta \f$ denotes the robots heading.
   */
  state_t state_estimate{state_t::Zero()};

  /**
   * @brief Recieved range and bearing measurement: 2x1 matrix.
   */
  measurement_t measurement{measurement_t::Zero()};

  /**
   * @brief The difference between the measurement recieved by the sensor and
   * the predicted measurement based on the vehicles' states.
   * @details The measurement residual is defined by the expression:
   * \f[ \tilde{\mathbf{y}}_k = \mathbf{z}_k - h(\hat{\mathbf{x}}_{k\mid
   * k-1})\f], where \f$\mathbf{z}\f$ is the measurement taken, and
   * \f$x_{k\mid k-1}\f$ is the estimated state of the robot.
   */
  measurement_t innovation{measurement_t::Zero()};

  /**
   * @brief Measurement innovation covariance matrix: 2x2 matrix.
   */
  measurementCovariance_t innovation_covariance{
      measurementCovariance_t::Zero()};

  /**
   * @brief Estimation Error Covariance: 3x3 matrix.
   * @details There is a high certainty in the prior value of system state,
   * therefore the prior estimation error covariance is initialised to a small
   * value.
   */
  covariance_t error_covariance{covariance_t::Identity() * 1e-3};

  /**
   * @brief Measurement noise covariance matrix: 2x2 matrix.
   * @details The matrix for the measurement noise covariance matrix take the
   * form \f[ v = \begin{bmatrix} q_r & 0 \\ 0 & q_\phi \end{bmatrix}, \f]
   * where \f$q_r\f$ denotes the range noise variance; and \f$\phi_r\f$
   * denotes the bearing noise variance.
   * @note The measurement noise is assumed to be uncorrelated and therefore
   * the covariance between the range and bearing is assumed to be zero.
   */
  measurementCovariance_t measurement_noise{measurementCovariance_t::Zero()};

  /**
   * @brief Information vector: 3x1 matrix.
   * @details Used by the information form of the (Extended) Kalman Filter.
   * The expression of the information vector take the form:
   * \f[\begin{align}\nabla &= \Sigma^{-1}\mathbf{x} \\ &= \Lambda
   * \mathbf{x}\end{align}\f], where \f$\Sigma\f$ denotes the estimation error
   * covariance matrix (Filter::EstimationParameters.error_covariance);
   * \f$\Lambda\f$ denotes the information matrix
   * (Filter::EstimationParameters.information_vector)
   * \f$\mathbf{x}\f$ denotes the state of the system
   * (Filter::EstimationParameters.state_estimate).
   */
  state_t information_vector{state_t::Zero()};

  /**
   * @brief Precision matrix: 3x3 matrix.
   * @details Used by the information form of the (Extended) Kalman Filter.
   * The expression for the information matrix takes the form \f[\Lambda =
   * \Sigma^{-1}\f], where \f$\Sigma\f$ denotes the estimation error
   * covariance matrix (Filter::EstimationParameters.error_covariance).
   */
  precision_t precision_matrix{error_covariance.inverse()};
};
} // namespace CL
