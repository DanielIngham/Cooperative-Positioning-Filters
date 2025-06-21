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

private:
  Eigen::Matrix<double, total_measurements, 1> accumulative_innovation =
      Eigen::Matrix<double, total_measurements, 1>::Zero();

  Eigen::Matrix<double, total_measurements, 1>
      accumulative_innovation_covariance =
          Eigen::Matrix<double, total_measurements, 1>::Zero();

  Eigen::Matrix<double, 2 + total_states, 1> accumulative_estimation_error =
      Eigen::Matrix<double, 2 + total_states, 1>::Zero();

  Eigen::Matrix<double, 2 + total_states, 1>
      accummulative_estimation_covariance =
          Eigen::Matrix<double, 2 + total_states, 1>::Zero();

  void prediction(const Robot::Odometry &, EstimationParameters &) override;

  void correction(EstimationParameters &, const EstimationParameters &,
                  const bool) override;

  void robustCorrection(EstimationParameters &, const EstimationParameters &);

  Eigen::Matrix<double, total_measurements, 1>
  computeMeasurementTau(const EstimationParameters &);

  Eigen::Matrix<double, 2 + total_states, 1>
  computeStateTau(const EstimationParameters &);
};

#endif // INCLUDE_INCLUDE_EKF_H_
