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

private:
  const unsigned short max_iterations = 2;

  void prediction(const Robot::Odometry &, EstimationParameters &) override;

  void correction(EstimationParameters &, const EstimationParameters &,
                  const bool) override;

  void robustCorrection(EstimationParameters &, const EstimationParameters &);
};
#endif // INCLUDE_INCLUDE_IEKF_H_
