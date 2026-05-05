#include "CLFilters/models/process.hpp"

namespace Filters::Models {
/**
 * @brief The unicycle motion model used to perform motion predictions.
 *
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The motion model used for robot takes the form:
 * \f[\begin{bmatrix} x_i^{(t+1)} \\  y_i^{(t+1)}
 * \\ \theta_i^{(t+1)}\end{bmatrix} = \begin{bmatrix} x_i^{(t)} +
 * \tilde{v}_i^{(t)}\Delta t\cos(\theta_i^{(t)}) \\ y_i^{(t)} +
 * \tilde{v}_i^{(t)}\Delta t\sin(\theta_i^{(t)})\\ \theta_i^{(t)} +
 * \tilde{\omega}_i^{(t)}\Delta t. \end{bmatrix}, \f] where \f$i\f$ denotes the
 * robots ID; \f$t\f$ denotes the current timestep; \f$\Delta t\f$ denotes the
 * sample period; \f$x\f$ and \f$y\f$ are the robots coordinates; \f$\theta\f$
 * denotes the robots heading (orientation); \f$\tilde{v}_i\f$ denotes the
 * forward velocity input; and \f$\tilde{\omega}\f$ denotes the angular velocity
 * input. Both \f$\tilde{v}_i\f$ and \f$\tilde{\omega}_i\f$ are normally
 * distributed random variables \f$\mathcal{N}(0,w)\f$ (see
 * Filter::EstimationParameters::process_noise).
 */
void Process::motionModel(const Data::Robot::Odometry &odometry, state_t &state,
                          const double sample_period) {

  state << state(X) + odometry.forward_velocity * sample_period *
                          std::cos(state(ORIENTATION)),
      state(Y) + odometry.forward_velocity * sample_period *
                     std::sin(state(ORIENTATION)),
      state(ORIENTATION) + odometry.angular_velocity * sample_period;
}

/**
 * @brief Calculates the Jacobian matrix of the unicycle motion model in terms
 * of the systems states (x,y,orientation).
 *
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The formula used for the calculation of the motion model
 * Jacobian takes the form:
 * \f[ F = \begin{bmatrix} 1 & 0 & -\tilde{v}\Delta t \sin(\theta) \\ 0 & 1
 * & \tilde{v} \Delta t \cos(\theta) \\ 0 & 0 & 1 \end{bmatrix}, \f]
 * where \f$\theta\f$ denotes the heading (orientation) of the ego vehicle;
 * and \f$\tilde{v}\f$ denotes the forward velocity. The forward velocity is
 * a random variable with Gaussian distributed noise \f$\mathcal{N}(0,w)\f$
 * , where \f$w\f$ is defined by the covariance matrix
 * Filter::EstimationParameters.measurement_noise. See Filter::motionModel for
 * information on the motion model from which this was derived.
 */
void Process::calculateMotionJacobian(
    const Data::Robot::Odometry &odometry,
    EstimationParameters &estimation_parameters, const double sample_period) {

  estimation_parameters.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, 1;
}

/**
 * @brief Calculates the Jacobian matrix of the motion model evaluated in terms
 * of the process inputs.
 *
 * @param[in,out] estimation_parameters The estimation parameters of the ego
 * robot.
 * @param[in] sample_period The period between odometry measurements.
 *
 * @details The formula used for the calculation of the process noise
 * Jacobian takes the form
 * \f[L = \begin{bmatrix}\Delta t \cos(\theta) & 0 \\ \Delta t \sin(\theta)
 * & 0 \\ 0 & \Delta t \end{bmatrix}, \f] where \f$\Delta t\f$ denotes the
 * sample period; and \f$\theta\f$ denotes the heading (orientation) of the
 * ego robot. measurement_noise. See EKF::prediction for information on the
 * motion model from which this was derived.
 */
void Process::calculateProcessJacobian(
    EstimationParameters &estimation_parameters, const double sample_period) {

  estimation_parameters.process_jacobian
      << sample_period *
             std::cos(estimation_parameters.state_estimate(ORIENTATION)),
      0,
      sample_period *
          std::sin(estimation_parameters.state_estimate(ORIENTATION)),
      0, 0, sample_period;
}

} // namespace Filters::Models
