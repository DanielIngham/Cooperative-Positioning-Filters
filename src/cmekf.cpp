#include "cmekf.h"
#include "Plotter.h"
#include "types.h"
#include <cmath>

namespace Filters {
CMEKF::CMEKF(Data::Handler &data) : EKF(data) {}

void CMEKF::correction(EstimationParameters &ego,
                       const EstimationParameters &agent) {

  const measurementJacobian_t agent_measurement_Jacobian{
      agentMeasurementJacobian(ego, agent)};

  measurementCovariance_t joint_sensor_noise{
      ego.measurement_noise + agent_measurement_Jacobian *
                                  agent.error_covariance *
                                  agent_measurement_Jacobian.transpose()};

  /* Covert the range and bearing into relative x and y coordinates.*/
  measurement_t measurement{relativePosition(ego.measurement)};

  /* Jacobian of the measurement with respect to the relative position change.*/
  const measurementCovariance_t J{jacobian(ego.measurement)};

  const measurement_t predicted_measurment{
      relativePosition(ego.state_estimate, agent.state_estimate)};

  /* Jacobian of the ego state respect to the relative position change.*/
  const measurementJacobian_t H{
      jacobian(ego.state_estimate, agent.state_estimate)};

  /* Calculate innovation Covariance. */
  ego.innovation_covariance = H * ego.error_covariance * H.transpose() +
                              J * joint_sensor_noise * J.transpose();

  /* Calculate Kalman Gain. */
  kalmanGain_t kalman_gain{ego.error_covariance * H.transpose() *
                           ego.innovation_covariance.inverse()};

  /* Calculate the innovation. */

  ego.innovation = measurement - predicted_measurment;
  Data::Robot::normaliseAngle(ego.innovation(BEARING));

  /* Update the state estimate. */
  ego.state_estimate += kalman_gain * ego.innovation;

  /* Update the estimation error covariance.  */
  ego.error_covariance -=
      kalman_gain * ego.innovation_covariance * kalman_gain.transpose();
}

Eigen::Matrix2d CMEKF::jacobian(Filters::measurement_t measurement) {
  measurementCovariance_t J;
  J(0, 0) = std::cos(measurement(BEARING));
  J(0, 1) = -measurement(RANGE) * std::sin(measurement(BEARING));

  J(1, 0) = std::sin(measurement(BEARING));
  J(1, 1) = measurement(RANGE) * std::cos(measurement(BEARING));
  return J;
}

Eigen::Matrix<double, 2, 3> CMEKF::jacobian(state_t ego, state_t agent) {
  const double delta_x{agent(X) - ego(X)};
  const double delta_y{agent(Y) - ego(Y)};

  Eigen::Matrix<double, 2, 3> H;
  H(0, 0) = -std::cos(ego(ORIENTATION));
  H(0, 1) = -std::sin(ego(ORIENTATION));
  H(0, 2) = -delta_x * std::sin(ego(ORIENTATION)) +
            delta_y * std::cos(ego(ORIENTATION));

  H(1, 0) = std::sin(ego(ORIENTATION));
  H(1, 1) = -std::cos(ego(ORIENTATION));
  H(1, 2) = -delta_x * std::cos(ego(ORIENTATION)) -
            delta_y * std::sin(ego(ORIENTATION));

  return H;
}

Eigen::Vector2d CMEKF::relativePosition(Filters::measurement_t measurement) {
  Eigen::Vector2d delta_local_pos;
  delta_local_pos(X) = measurement(RANGE) * std::cos(measurement(BEARING));
  delta_local_pos(Y) = measurement(RANGE) * std::sin(measurement(BEARING));
  return delta_local_pos;
}

Eigen::Vector2d CMEKF::relativePosition(state_t ego, state_t agent) {

  Eigen::Vector2d delta_global_pos{agent.head(2) - ego.head(2)};
  Eigen::Rotation2Dd Rot(-ego(ORIENTATION));

  return Rot * delta_global_pos;
}

} // namespace Filters
