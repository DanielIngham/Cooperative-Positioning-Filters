/**
 * @file filter.h
 * @brief Header file of the parent class containing shared functionality amoung
 * cooperative localsiation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/sensors/odom_data.hpp"

#include <Eigen/Dense>

namespace CL::filter {
class Filter {
public:
  Filter() = delete;
  Filter(Filter &&) = default;
  Filter(const Filter &) = default;
  Filter &operator=(Filter &&) = default;
  Filter &operator=(const Filter &) = default;
  virtual ~Filter() = default;

  explicit Filter(const EstimationParameters &prior) : prior_{prior} {};

  [[nodiscard]] virtual EstimationParameters
  prediction(sensors::OdomData const &odometry, EstimationParameters const &ego,
             double sample_period) = 0;

  virtual void correction(EstimationParameters &ego,
                          const EstimationParameters &agent) = 0;

private:
protected:
  EstimationParameters prior_;
};
} // namespace CL::filter
