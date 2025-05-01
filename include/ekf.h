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
  static const unsigned short total_inputs = 3;

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
    double state_estimate[total_states][1] = {{0.0}, {0.0}, {0.0}};

    /**
     * @brief Estimation Error Covariance.
     */
    double error_covarince[total_states][total_states] = {
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    double kalman_gain[total_states][total_states] = {
        {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    /**
     * @brief Odometry process noise: forward velocity, angular velocity.
     * @details \f[v & \omega \f]
     */
    double process_noise[total_inputs][1] = {{0.0}, {0.0}};

    /**
     * @brief Measurement noise: range, bearing.
     * @details \f[r & \phi \f]
     */
    double measurement_noise[total_inputs][1] = {{0.0}, {0.0}};

    double motion_jacobian[total_states][total_states];

    double measurment_jacobian[total_states][total_states];
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
