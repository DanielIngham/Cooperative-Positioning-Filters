#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

namespace CL::Models {
class Measurement {
public:
  Measurement() = default;
  Measurement(Measurement &&) = default;
  Measurement(const Measurement &) = default;
  Measurement &operator=(Measurement &&) = default;
  Measurement &operator=(const Measurement &) = default;
  ~Measurement() = default;

  Measurement(const EstimationParameters &ego,
              const EstimationParameters agent);

  [[nodiscard]] static measurement_t measurementModel(const state_t &,
                                                      const state_t &);

  [[nodiscard]] static double rangeMeasurementModel(const state_t &,
                                                    const state_t &);

  void calculateMeasurementJacobian(const EstimationParameters &,
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
  /**
   * Getter for the ego vehicle's measurement Jacobian matrix.
   * @returns a Jacobian matrix.
   */
  measurementJacobian_t getEgoJacobian() const;

  /**
   * Getter for the observed agents measurement Jacobian matrix.
   * @returns a Jacobian matrix;
   */
  measurementJacobian_t getAgentJacobian() const;

  /**
   * Getter for the Jacobian matrix that contains both the ego vehicle and agent
   * Jacobian matrices.
   */
  augmentedMeasurementJacobian_t getAugmentedJacobian() const;

  /**
   * Getter for the predicted measurement given the state of the ego vehicle and
   * the observed agent.
   */
  measurement_t getPrediction() const;

private:
  measurement_t predicted_measurement_{};

  /**
   * @brief Jacobian of the measurement model: 2 x 6 matrix.
   */
  augmentedMeasurementJacobian_t measurement_jacobian_{};
};
} // namespace CL::Models
