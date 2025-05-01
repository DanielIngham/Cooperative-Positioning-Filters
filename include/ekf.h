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
 */
class EKF {
public:
  EKF(DataHandler &data);
  ~EKF();
  void peformInference();

private:
  static const unsigned short total_states = 3;
  static const unsigned short total_inputs = 2;
  static const unsigned short total_measurements = 2;

  enum { FORWARD_VELOCITY = 0, ANGULAR_VELOCITY = 1 };
  enum { RANGE = 0, BEARING = 1 };
  enum { X = 0, Y = 1, ORIENTATION = 2 };

  DataHandler &data_;

  /**
   * @struct EstimationParameters
   * @brief Houses the parameters required for performing Bayesian filtering.
   */
  struct EstimationParameters {

    /**
     * @brief Estimated robot state: x-coordinate, y-coordinate, orientation.
     * @details \f[begin{bmatrix} x & y & \theta
     * \end{bmatrix}^\top\f]
     */
    Eigen::Matrix<double, total_states, 1> state_estimate =
        Eigen::Matrix<double, total_states, 1>::Zero();

    /**
     * @brief Estimation Error Covariance.
     */
    Eigen::Matrix<double, total_states, total_states> error_covarince =
        Eigen::Matrix<double, total_states, total_states>::Zero();

    /**
     * @brief Kalman gain.
     */
    Eigen::Matrix<double, total_states, total_states> kalman_gain =
        Eigen::Matrix<double, total_states, total_states>::Zero();

    /**
     * @brief Odometry process noise: forward velocity, angular velocity.
     * @details \f[v & \omega \f]
     */
    Eigen::Matrix<double, total_inputs, 1> process_noise =
        Eigen::Matrix<double, total_inputs, 1>::Zero();

    /**
     * @brief Measurement noise: range, bearing.
     * @details \f[r & \phi \f]
     */
    Eigen::Matrix<double, total_measurements, 1> measurement_noise =
        Eigen::Matrix<double, total_measurements, 1>::Zero();

    /**
     * @brief Jacobian of the motion model.
     */
    Eigen::Matrix<double, total_states, total_states> motion_jacobian =
        Eigen::Matrix<double, total_states, total_states>::Zero();

    /**
     * @brief Jacobian of the measurement model.
     */
    Eigen::Matrix<double, total_measurements, total_measurements>
        measurment_jacobian = Eigen::Matrix<double, total_measurements,
                                            total_measurements>::Zero();
  };

  void prediction(const Robot::Odometry &odometry,
                  const Robot::State &prior_state,
                  EstimationParameters &estimation_parameters);

  void correction(const Robot::Odometry &odometry,
                  const Robot::State prior_state,
                  const EstimationParameters estimation_parameters);
  /**
   * @brief Houses all estimation parameters for all robots.
   */
  std::vector<EstimationParameters> robots;
};

#endif // INCLUDE_INCLUDE_EKF_H_
