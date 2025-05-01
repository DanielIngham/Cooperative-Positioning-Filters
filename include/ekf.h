#ifndef INCLUDE_INCLUDE_EKF_H_
#define INCLUDE_INCLUDE_EKF_H_

#include "data_handler.h"
#include <vector>

/**
 * @class EKF
 * @brief Extended Kalman Filter.
 */
class EKF {
public:
  EKF(DataHandler &data);
  ~EKF();

  void setDataSet();
  void peformInference();

private:
  DataHandler &data_;

  void prediction();
  void correction();

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
    double state_estimate[3][1];

    /**
     * @brief Estimation Error Covariance.
     */
    double error_covarince[3][3];

    double kalman_gain[3][3];

    /**
     * @brief Odometry process noise: forward velocity, angular velocity.
     * @details \f[v & \omega \f]
     */
    double process_noise[2][1];

    /**
     * @brief Measurement noise: range, bearing.
     * @details \f[r & \phi \f]
     */
    double measurement_noise[2][1];
  };

  /*
   * @brief Houses all estimation parameters for all robots.
   */
  std::vector<EstimationParameters> robots;
};

#endif // INCLUDE_INCLUDE_EKF_H_
