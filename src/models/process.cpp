#include "CL/models/process.hpp"
#include "CL/common/types.hpp"
#include "CL/utils/utils.hpp"

namespace CL::Models {
Process::Process(sensors::OdomData const &odometry, state_t const &state,
                 const double sample_period) {

  predicted_state_ = motionModel(odometry, state, sample_period);
  motion_jacobian_ = calculateMotionJacobian(odometry, state, sample_period);
  process_jacobian_ = calculateProcessJacobian(state, sample_period);
}

state_t Process::predictedState() { return predicted_state_; }

motionJacobian_t Process::motionJacobian() { return motion_jacobian_; }

processJacobian_t Process::processJacobian() { return process_jacobian_; }

state_t Process::motionModel(sensors::OdomData const &odometry,
                             state_t const &state, double const sample_period,
                             input_t noise) {

  state_t predicted_state{state};
  input_t const input{odometry.input() + noise};

  predicted_state(X) = state(X) + input(FORWARD_VELOCITY) * sample_period *
                                      std::cos(state(ORIENTATION));
  predicted_state(Y) = state(Y) + input(FORWARD_VELOCITY) * sample_period *
                                      std::sin(state(ORIENTATION));
  predicted_state(ORIENTATION) =
      state(ORIENTATION) + input(ANGULAR_VELOCITY) * sample_period;

  utils::normaliseAngle(predicted_state(ORIENTATION));

  return predicted_state;
}

motionJacobian_t
Process::calculateMotionJacobian(sensors::OdomData const &odometry,
                                 state_t const &state,
                                 const double sample_period) {

  motionJacobian_t motion_jacobian;
  motion_jacobian(X, X) = 1;
  motion_jacobian(X, Y) = 0;
  motion_jacobian(X, ORIENTATION) = -odometry.input(FORWARD_VELOCITY) *
                                    sample_period *
                                    std::sin(state(ORIENTATION));

  motion_jacobian(Y, X) = 0;
  motion_jacobian(Y, Y) = 1;
  motion_jacobian(Y, ORIENTATION) = odometry.input(FORWARD_VELOCITY) *
                                    sample_period *
                                    std::cos(state(ORIENTATION));

  motion_jacobian(ORIENTATION, X) = 0;
  motion_jacobian(ORIENTATION, Y) = 0;
  motion_jacobian(ORIENTATION, ORIENTATION) = 1;

  return motion_jacobian;
}

processJacobian_t
Process::calculateProcessJacobian(const state_t &state,
                                  const double sample_period) {

  processJacobian_t process_jacobian{};
  process_jacobian(X, FORWARD_VELOCITY) =
      sample_period * std::cos(state(ORIENTATION));
  process_jacobian(X, ANGULAR_VELOCITY) = 0;

  process_jacobian(Y, FORWARD_VELOCITY) =
      sample_period * std::sin(state(ORIENTATION));
  process_jacobian(Y, ANGULAR_VELOCITY) = 0;

  process_jacobian(ORIENTATION, FORWARD_VELOCITY) = 0;
  process_jacobian(ORIENTATION, ANGULAR_VELOCITY) = sample_period;

  return process_jacobian;
}

} // namespace CL::Models
