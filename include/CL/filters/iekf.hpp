#ifndef INCLUDE_INCLUDE_IEKF_H_
#define INCLUDE_INCLUDE_IEKF_H_

/**
 * @file iekf.h
 * @brief Header file of the Iterative Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-20
 */
#include "CL/filters/ekf.hpp"

#include <Eigen/Dense>
#include <UtiasMrclam/DataHandler.hpp>

namespace CL::filter {
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
  IEKF() = delete;
  IEKF(IEKF &&) = default;
  IEKF(const IEKF &) = default;
  IEKF &operator=(IEKF &&) = delete;
  IEKF &operator=(const IEKF &) = delete;
  ~IEKF() = default;

  IEKF(const EstimationParameters &prior) : EKF{prior} {};

  /**
   * @brief Performs the Iterative Extended Kalman correct step.
   * @param[in,out] ego_robot The estimation parameters of the ego robot.
   * @param[in] other_agent The estimation parameters of the obejct that was
   * measured by the ego robot.
   */
  void correction(EstimationParameters &ego, const EstimationParameters &agent,
                  sensors::MeasData const &meas) override;

  /**
   * @brief A robust version of the correction function that uses the Huber cost
   * function to increase estimation error covariance of measurements that seem
   * to be outliers.
   * @param[in,out] ego_robot The estimation parameters of the ego robot.
   * @param[in] other_agent The estimation parameters of the obejct that was
   * measured by the ego robot.
   */
  void robustCorrection(EstimationParameters &, const EstimationParameters &);

private:
  /** Maximum number of iterations performed by the iterative update. */
  const unsigned short max_iterations_{100U};
};
} // namespace CL::filter
#endif // INCLUDE_INCLUDE_IEKF_H_
