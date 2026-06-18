/**
 * @file information_filter.h
 * @brief Header file of the Information form for the Extended Kalman Fitler
 * implementation for multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#ifndef INCLUDE_INCLUDE_INFORMATION_H_
#define INCLUDE_INCLUDE_INFORMATION_H_

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include "CL/filters/filter.hpp"

#include <UtiasMrclam/DataHandler.hpp>

namespace CL::filter {

/**
 * @class InformationFilter
 * Information form of the Extended Kalman Filter.
 */
class InformationFilter : public Filter {
public:
  /**
   * @brief InformationFilter class constructor.
   * @details This constructor sets up the prior states and parameters to
   * perform Extended Information filtering.
   * @param[in] data Class containing all robot and landmark data.
   */
  InformationFilter() = delete;
  InformationFilter(InformationFilter &&) = default;
  InformationFilter(const InformationFilter &) = default;
  InformationFilter &operator=(InformationFilter &&) = delete;
  InformationFilter &operator=(const InformationFilter &) = delete;
  ~InformationFilter() = default;

  InformationFilter(const EstimationParameters &prior) : Filter{prior} {};

  /**
   * @brief performs the prediction step of the Information filter.
   * @param[in] odometry The prior inputs into the system comprising a forward
   * and angular velocity.
   * @param[in,out] parameters The parameters required by the
   * Information filter to perform the prediction step.
   */
  EstimationParameters prediction(sensors::OdomData const &odometry,
                                  EstimationParameters const &parameters,
                                  double sample_period) override;

  /**
   * @brief Performs Information Filter correct step.
   * @param[in,out] ego_robot The parameters required by the Extended
   * Kalman filter to perform the correction step.
   * @param[in] other_agent The robot that was measured by the ego robot.
   * @param[in] robust Flag which determines whether the information and
   * precision should be updated using a robust cost function.
   */
  void correction(EstimationParameters &ego, EstimationParameters const &agent,
                  sensors::MeasData const &meas) override;

private:
};
} // namespace CL::filter
#endif // INCLUDE_INCLUDE_INFORMATION_H_
