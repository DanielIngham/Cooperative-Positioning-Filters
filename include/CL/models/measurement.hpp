#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

namespace CL::Models {
class Measurement {
public:
  Measurement() = default;
  Measurement(Measurement &&) = delete;
  Measurement(const Measurement &) = delete;
  Measurement &operator=(Measurement &&) = delete;
  Measurement &operator=(const Measurement &) = delete;
  ~Measurement() = default;

  [[nodiscard]] static measurement_t measurementModel(const state_t &,
                                                      const state_t &);

  [[nodiscard]] static double rangeMeasurementModel(const state_t &,
                                                    const state_t &);

  static void calculateMeasurementJacobian(EstimationParameters &,
                                           const EstimationParameters &);

  [[nodiscard]] static measurementJacobian_t
  egoMeasurementJacobian(const EstimationParameters &,
                         const EstimationParameters &);

  [[nodiscard]] static measurementJacobian_t
  agentMeasurementJacobian(const EstimationParameters &,
                           const EstimationParameters &);

  [[nodiscard]] static vector3D_t
  egoRangeMeasurementJacobian(const EstimationParameters &,
                              const EstimationParameters &);

  [[nodiscard]] static vector3D_t
  agentRangeMeasurementJacobian(const EstimationParameters &,
                                const EstimationParameters &);

private:
};
} // namespace CL::Models
