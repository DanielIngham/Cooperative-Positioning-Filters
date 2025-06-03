#ifndef INCLUDE_INCLUDE_IEKF_H_
#define INCLUDE_INCLUDE_IEKF_H_

/**
 * @file iekf.h
 * @brief Header file of the Iterative Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-20
 */
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
class IEKF : public Filter {
public:
  explicit IEKF(DataHandler &data);
  ~IEKF() override;
  void performInference();

private:
  void prediction(const Robot::Odometry &, EstimationParameters &);

  void correction(EstimationParameters &, const EstimationParameters &);
  void robustCorrection(EstimationParameters &, const EstimationParameters &);

  Eigen::Matrix<double, total_measurements, total_measurements>
  HuberMeasurement(const Eigen::Matrix<double, total_measurements, 1> &);

  Eigen::Matrix<double, total_states + 2, total_states + 2>
  HuberState(const Eigen::Matrix<double, total_states + 2, 1> &);
};
#endif // INCLUDE_INCLUDE_IEKF_H_
