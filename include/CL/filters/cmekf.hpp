#pragma once
#include "CL/common/types.hpp"
#include "CL/filters/ekf.hpp"

#include <Eigen/src/Core/Matrix.h>
#include <UtiasMrclam/DataHandler.hpp>

namespace CL::filter {
class CMEKF : public EKF {
public:
  CMEKF() = delete;
  CMEKF(CMEKF &&) = default;
  CMEKF(const CMEKF &) = default;
  CMEKF &operator=(CMEKF &&) = delete;
  CMEKF &operator=(const CMEKF &) = delete;
  ~CMEKF() = default;

  CMEKF(const EstimationParameters &prior) : EKF{prior} {};

  void correction(EstimationParameters &ego, EstimationParameters const &agent,
                  sensors::MeasData const &meas) override;

protected:
  Eigen::Matrix2d jacobian(measurement_t measurement);
  Eigen::Matrix<double, total_measurements, total_states>
  jacobian(state_t ego, state_t agent);

  Eigen::Vector2d relativePosition(measurement_t measurement);
  Eigen::Vector2d relativePosition(state_t ego, state_t agent);
};

} // namespace CL::filter
