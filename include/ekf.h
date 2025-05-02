/**
 * @file ekf.h
 * @brief Header file of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#ifndef INCLUDE_INCLUDE_EKF_H_
#define INCLUDE_INCLUDE_EKF_H_

#include <DataHandler/DataHandler.h>
#include <Eigen/Dense>
#include <vector>

/**
 * @class EKF
 * @brief Extended Kalman Filter.
 * @details The data contained in the class includes:
 * - All robots states, odometry readings, and measurements.
 * - All landmarks positions.
 * - All robot sensors errors statistics.
 */
class EKF {
public:
  EKF(DataHandler &data);
  ~EKF();
  void peformInference();

private:
  /**
   * @brief The total number of system states: x coordinate [m], y-coordinate
   * [m], and orientation (heading) [rad].
   */
  static const unsigned short total_states = 3;
  /**
   * @brief The total number of system inputs: forward velocity [m/s] and
   * angular velocity [rad/s].
   */
  static const unsigned short total_inputs = 2;
  /**
   * @brief The total number of system measurements: range [m] and bearing
   * [rad].
   */
  static const unsigned short total_measurements = 2;

  enum { FORWARD_VELOCITY = 0, ANGULAR_VELOCITY = 1 };
  enum { RANGE = 0, BEARING = 1 };
  enum { X = 0, Y = 1, ORIENTATION = 2 };

  /**
   * @brief Data class housing all the data pertaining to cooperative
   * localisation (positioning).
   */
  DataHandler &data_;

  /**
   * @struct EstimationParameters
   * @brief Houses the parameters required for performing Bayesian filtering.
   */
  struct EstimationParameters {

    /**
     * @brief Estimated robot state - x coordinate [m], y-coordinate [m], and
     * orientation (heading) [rad]: 3x1 matrix.
     * @details The state vector of the robot take the form \f[\begin{bmatrix} x
     * & y & \theta \end{bmatrix}^\top, \f] where \f$x\f$ and \f$y\f$ denotes
     * the robots 2D coordinates; and \f$\theta \f$ denotes the robots heading.
     */
    Eigen::Matrix<double, total_states, 1> state_estimate =
        Eigen::Matrix<double, total_states, 1>::Zero();

    /**
     * @brief Recieved range and bearing measurement: 2x1 matrix.
     */
    Eigen::Matrix<double, total_measurements, 1> measurement =
        Eigen::Matrix<double, total_measurements, 1>::Zero();

    /**
     * @brief Estimation Error Covariance: 3x3 matrix.
     * @details There is a high certainty in the prior value of system state,
     * therefore the prior estimation error covariance is initialised to a small
     * value.
     */
    Eigen::Matrix<double, total_states, total_states> error_covariance =
        Eigen::Matrix<double, total_states, total_states>::Identity() * 0.001;

    /**
     * @brief Measurement correction innovation matrix: 2x2 matrix.
     */
    Eigen::Matrix<double, total_measurements, total_measurements> innovation =
        Eigen::Matrix<double, total_measurements, total_measurements>::Zero();

    /**
     * @brief Kalman gain: 6x2 matrix.
     * @note The reason the Kalman gain matrix has 6 elements is because the
     * cooperative localisation (positioning) requires the estimation
     * covariances of both the ego and measured robots [3x2]. See
     * EKF::correction.
     */
    Eigen::Matrix<double, 2 * total_states, total_measurements> kalman_gain =
        Eigen::Matrix<double, 2 * total_states, total_measurements>::Zero();

    /**
     * @brief Odometry process noise covariance matrix: 2x2 matrix.
     * @details The process noise covariance matrix is defined by the
     * expression:
     * \f[ w = \begin{bmatrix} q_v & 0 \\ 0 & q_\omega \end{bmatrix}, \f] where
     * \f$q_v\f$ denotes the forward velocity noise variance; and
     * \f$q_\omega\f$ denotes the angular velocity noise variance.
     * @note The process noise is assumed to be uncorrelated and therefore the
     * covariance between the forward velocity and the angular velocity is
     * assumed to be zero.
     */
    Eigen::Matrix<double, total_inputs, total_inputs> process_noise =
        Eigen::Matrix<double, total_inputs, total_inputs>::Zero();

    /**
     * @brief Measurement noise covariance matrix: 2x2 matrix.
     * @details The matrix for the measurement noise covariance matrix take the
     * form \f[ v = \begin{bmatrix} q_r & 0 \\ 0 & q_\phi \end{bmatrix}, \f]
     * where \f$q_r\f$ denotes the range noise variance; and \f$\phi_r\f$
     * denotes the bearing noise variance.
     * @note The measurement noise is assumed to be uncorrelated and therefore
     * the covariance between the range and bearing is assumed to be zero.
     */
    Eigen::Matrix<double, total_measurements, total_measurements>
        measurement_noise = Eigen::Matrix<double, total_measurements,
                                          total_measurements>::Zero();

    /**
     * @brief Jacobian matrix of the motion model.
     * @details The formula used for the calculation of the motion model
     * Jacobian takes the form:
     * \f[ F = \begin{bmatrix} 1 & 0 & -\tilde{v}\Delta t \sin(\theta) \\ 0 & 1
     * & \tilde{v} \Delta t \cos(\theta) \\ 0 & 0 & 1 \end{bmatrix}, \f]
     * where \f$\theta\f$ denotes the heading (orientation) of the ego vehicle;
     * and \f$\tilde{v}\f$ denotes the forward velocity. The forward velocity is
     * a random variable with Gaussian distributed noise \f$\mathcal{N}(0,w)\f$
     * , where \f$w\f$ is defined by the covariance matrix
     * EKF::EstimationParameters.measurement_noise. See EKF::prediction for
     * information on the motion model from which this was derived.
     */
    Eigen::Matrix<double, total_states, total_states> motion_jacobian =
        Eigen::Matrix<double, total_states, total_states>::Zero();

    /**
     * @brief Jacobian matrix of the process noise.
     * @details The formula used for the calculation of the process noise
     * Jacobian takes the form
     * \f[L = \begin{bmatrix}\Delta t \cos(\theta) & 0 \\ \Delta t \sin(\theta)
     * & 0 \\ 0 & \Delta t \end{bmatrix}, \f] where \f$\Delta t\f$ denotes the
     * sample period; and \f$\theta\f$ denotes the heading (orientation) of the
     * ego robot. measurement_noise. See EKF::prediction for information on the
     * motion model from which this was derived.
     */
    Eigen::Matrix<double, total_states, total_inputs> process_jacobian =
        Eigen::Matrix<double, total_states, total_inputs>::Zero();
    /**
     * @brief Jacobian of the measurement model.
     * @details The formula used for the calculation of the Jacobian of the
     * measurement matrix between ego vehicle \f$i\f$ and measured vehicle
     * \f$j\f$ take the form
     * \f[ H = \begin{bmatrix} \frac{-\Delta x}{d} & \frac{-\Delta y}{d} & 0 &
     * \frac{\Delta x}{d} & \frac{\Delta y}{d} & 0 \\ \frac{\Delta y}{d^2} &
     * \frac{-\Delta x}{d^2} & -1 & \frac{-\Delta y}{d^2} & \frac{\Delta x}{d^2}
     * & 0 \end{bmatrix} \f] where \f$\Delta x = x_j - x_i\f$; \f$\Delta y = y_j
     * - y_i\f$; and \f$\Delta d = \sqrt{\Delta x^2 + \Delta y^2}\f$.
     */
    Eigen::Matrix<double, total_measurements, 2 * total_states>
        measurment_jacobian =
            Eigen::Matrix<double, total_measurements, 2 * total_states>::Zero();
  };

  void prediction(const Robot::Odometry &, EstimationParameters &);

  void correction(EstimationParameters &, const EstimationParameters &);
  /**
   * @brief Houses all estimation parameters for all robots.
   */
  std::vector<EstimationParameters> robot_parameters;
};

#endif // INCLUDE_INCLUDE_EKF_H_
