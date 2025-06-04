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
  unsigned long total_observations = 0;

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

  void prediction(const Robot::Odometry &, EstimationParameters &);

  void correction(EstimationParameters &, const EstimationParameters &,
                  const bool);

  void robustCorrection(EstimationParameters &, const EstimationParameters &,
                        Eigen::Matrix<double, total_measurements, 1> &,
                        Eigen::Matrix<double, 2 + total_states, 1> &);

  Eigen::Matrix<double, total_measurements, 1>
  computeMeasurementTau(const EstimationParameters &);

  Eigen::Matrix<double, 2 + total_states, 1>
  computeStateTau(const EstimationParameters &);

  Eigen::Matrix<double, total_measurements, total_measurements>
  HuberMeasurement(const Eigen::Matrix<double, total_measurements, 1> &,
                   const Eigen::Matrix<double, total_measurements, 1> &);

  Eigen::Matrix<double, 2 + total_states, total_states + 2>
  HuberState(const Eigen::Matrix<double, 2 + total_states, 1> &,
             const Eigen::Matrix<double, 2 + total_states, 1> &);

  Eigen::Matrix<double, total_measurements, 1>
  calculateNormalisedMeasurementResidual(
      const Eigen::Matrix<double, total_measurements, total_measurements> &,
      const Eigen::Matrix<double, total_measurements, 1> &);

  Eigen::Matrix<double, 2 + total_states, 1>
  calculateNormalisedEstimationResidual(
      const Eigen::Matrix<double, 2 + total_states, total_measurements> &,
      const Eigen::Matrix<double, total_measurements, total_measurements> &,
      const Eigen::Matrix<double, 2 + total_states, 1> &);
};

#endif // INCLUDE_INCLUDE_EKF_H_
