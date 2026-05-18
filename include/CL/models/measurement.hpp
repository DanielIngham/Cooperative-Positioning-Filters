#pragma once

#include "CL/common/estimation_parameters.hpp"
#include <type_traits>

#include <Eigen/Dense>

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
  Eigen::MatrixXd getEgoJacobian() const;

  /**
   * Getter for the observed agents measurement Jacobian matrix.
   * @returns a Jacobian matrix;
   */
  Eigen::MatrixXd getAgentJacobian() const;

  /**
   * Getter for the Jacobian matrix that contains both the ego vehicle and agent
   * Jacobian matrices.
   */
  Eigen::MatrixXd getAugmentedJacobian() const;

  /**
   * Getter for the predicted measurement given the state of the ego vehicle and
   * the observed agent.
   */
  Eigen::MatrixXd getPrediction() const;

protected:
  Measurement() = default;

  Eigen::MatrixXd predicted_measurement_{};

  /**
   * @brief Jacobian of the measurement model: 2 x 6 matrix.
   */
  Eigen::MatrixXd measurement_jacobian_{};

  /**
   * Jacobian matrix of the states of the ego robot with respect to the
   * measurement model.
   */
  Eigen::MatrixXd ego_jacobian_{};

  /**
   * Jacobian matrix of the states of the observed agent with respect to the
   * measurement model.
   */
  Eigen::MatrixXd agent_jacobian_{};

private:
};
} // namespace CL::Models
