/**
 * @file ekf.h
 * @brief Header file of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#ifndef INCLUDE_INCLUDE_EKF_H_
#define INCLUDE_INCLUDE_EKF_H_

#include "CL/common/estimation_parameters.hpp"
#include "CL/filters/filter.hpp"

#include <Eigen/Dense>
#include <UtiasMrclam/DataHandler.hpp>

namespace CL::filter {

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
  EKF() = delete;
  EKF(EKF &&) = default;
  EKF(const EKF &) = default;
  EKF &operator=(EKF &&) = delete;
  EKF &operator=(const EKF &) = delete;
  ~EKF() = default;

  EKF(const EstimationParameters &prior) : Filter{prior} {};

  [[nodiscard]] EstimationParameters
  prediction(const Data::Robot::Odometry &odometry,
             const EstimationParameters &parameters,
             double sample_period) override;

  void correction(EstimationParameters &ego,
                  const EstimationParameters &agent) override;

private:
  void robustCorrection(EstimationParameters &, const EstimationParameters &);
};
} // namespace CL::filter
#endif // INCLUDE_INCLUDE_EKF_H_
