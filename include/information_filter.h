/**
 * @file information_filter.h
 * @brief Header file of the Information form for the Extended Kalman Fitler
 * implementation for multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#ifndef INCLUDE_INCLUDE_INFORMATION_H_
#define INCLUDE_INCLUDE_INFORMATION_H_

#include "filter.h"

#include <DataHandler.h>

namespace Filter {

/**
 * @class InformationFilter
 * Information form of the Extended Kalman Filter.
 */
class InformationFilter : public Filter {
public:
  InformationFilter(Data::Handler &data);
  InformationFilter(InformationFilter &&) = default;
  InformationFilter(const InformationFilter &) = default;
  InformationFilter &operator=(InformationFilter &&) = delete;
  InformationFilter &operator=(const InformationFilter &) = delete;
  ~InformationFilter() override;

  void prediction(const Data::Robot::Odometry &,
                  EstimationParameters &) override;

  void correction(EstimationParameters &,
                  const EstimationParameters &) override;

private:
};
} // namespace Filter
#endif // INCLUDE_INCLUDE_INFORMATION_H_
