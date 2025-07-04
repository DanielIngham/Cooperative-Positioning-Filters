#ifndef INCLUDE_INCLUDE_IEKF_H_
#define INCLUDE_INCLUDE_IEKF_H_

/**
 * @file iekf.h
 * @brief Header file of the Iterative Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-20
 */
#include "ekf.h"

#include <DataHandler/DataHandler.h>
#include <Eigen/Dense>

/**
 * @class IEKF
 * @brief Itertative Extended Kalman Filter (IEKF).
 * @details This is a child of the EKF since it has the same prediction step.
 * The difference lies in the Gauss-Nwton optimisation approach to the
 * measurement update.
 * @note This filter does not work well in the precense of outliers and
 * therefore it is recommended to use the robustCorrection. Even so, since the
 * system does not exhibit strong non-linearities, the benifits of relinearising
 * the Jacobians is minimial. In conclusion, using the standard EKF (or EKF with
 * robust correction) results in better performance.
 */
class IEKF : public EKF {
public:
  explicit IEKF(DataHandler &data);
  ~IEKF() override;

private:
  const unsigned short max_iterations = 1;

  void correction(EstimationParameters &,
                  const EstimationParameters &) override;

  void robustCorrection(EstimationParameters &, const EstimationParameters &);
};
#endif // INCLUDE_INCLUDE_IEKF_H_
