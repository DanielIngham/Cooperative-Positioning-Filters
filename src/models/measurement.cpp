#include "CL/models/measurement.hpp"
#include <Eigen/src/Core/Matrix.h>

namespace CL::Models {

Eigen::MatrixXd Measurement::getEgoJacobian() const { return ego_jacobian_; }

Eigen::MatrixXd Measurement::getAgentJacobian() const {
  return agent_jacobian_;
}

Eigen::MatrixXd Measurement::getAugmentedJacobian() const {
  return measurement_jacobian_;
}

Eigen::MatrixXd Measurement::getPrediction() const {
  return predicted_measurement_;
}

} // namespace CL::Models
