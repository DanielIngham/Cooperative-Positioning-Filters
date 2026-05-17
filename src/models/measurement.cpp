#include "CL/models/measurement.hpp"

namespace CL::Models {

measurementJacobian_t Measurement::getEgoJacobian() const {
  return measurement_jacobian_.leftCols<3>();
}

measurementJacobian_t Measurement::getAgentJacobian() const {
  return measurement_jacobian_.rightCols<3>();
}

augmentedMeasurementJacobian_t Measurement::getAugmentedJacobian() const {
  return measurement_jacobian_;
}

measurement_t Measurement::getPrediction() const {
  return predicted_measurement_;
}

} // namespace CL::Models
