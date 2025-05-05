/**
 * @file ekf.h
 * @brief Header file of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#ifndef INCLUDE_INCLUDE_EKF_H_
#define INCLUDE_INCLUDE_EKF_H_

#include "filter.h"

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
class EKF : public Filter {
public:
  explicit EKF(DataHandler &data);
  ~EKF() override;
  void performInference();

private:
  void prediction(const Robot::Odometry &, EstimationParameters &);

  void correction(EstimationParameters &, const EstimationParameters &);
};

#endif // INCLUDE_INCLUDE_EKF_H_
