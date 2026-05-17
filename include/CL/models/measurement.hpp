#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include <type_traits>

namespace CL::Models {
class Measurement {
public:
  Measurement(Measurement &&) = default;
  Measurement(const Measurement &) = default;
  Measurement &operator=(Measurement &&) = default;
  Measurement &operator=(const Measurement &) = default;
  ~Measurement() = default;

  template <typename T>
  [[nodiscard]] static T generateMeasurement(const EstimationParameters &ego,
                                             const EstimationParameters agent) {
    static_assert(std::is_base_of_v<Measurement, T>,
                  "T must be child of Measurement class.");
    T measurement{ego, agent};
    return measurement;
  };

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

protected:
  Measurement() = default;

  measurement_t predicted_measurement_{};

  /**
   * @brief Jacobian of the measurement model: 2 x 6 matrix.
   */
  augmentedMeasurementJacobian_t measurement_jacobian_{};

  virtual measurement_t model(const state_t &, const state_t &) = 0;

private:
};
} // namespace CL::Models
