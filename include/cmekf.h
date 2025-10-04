#pragma once
#include "ekf.h"
#include "types.h"

#include <DataHandler.h>
#include <Eigen/src/Core/Matrix.h>

namespace Filters {
class CMEKF : public EKF {
public:
  CMEKF(Data::Handler &data);
  CMEKF(CMEKF &&) = default;
  CMEKF(const CMEKF &) = default;
  CMEKF &operator=(CMEKF &&) = delete;
  CMEKF &operator=(const CMEKF &) = delete;
  ~CMEKF() = default;

  void correction(EstimationParameters &,
                  const EstimationParameters &) override;

protected:
  Eigen::Matrix2d jacobian(measurement_t measurement);
  Eigen::Matrix<double, total_measurements, total_states>
  jacobian(state_t ego, state_t agent);

  Eigen::Vector2d relativePosition(measurement_t measurement);
  Eigen::Vector2d relativePosition(state_t ego, state_t agent);
};

} // namespace Filters
