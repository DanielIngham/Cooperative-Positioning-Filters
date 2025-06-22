/**
 * @file iekf.cpp
 * @brief Implementation of the Iterative Extended Kalman Fitler implementation
 * for multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#include "iekf.h"
#include "filter.h"

#include <DataHandler/Landmark.h>
#include <DataHandler/Robot.h>
#include <Eigen/Cholesky>
#include <stdexcept>

/**
 * @brief IEKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Iterative Extended Kalman filtering.
 * @param[in] data Class containing all robot and landmark data.
 */
IEKF::IEKF(DataHandler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
IEKF::~IEKF() {}

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 */
void IEKF::prediction(const Robot::Odometry &odometry,
                      EstimationParameters &estimation_parameters) {

  const double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, estimation_parameters, sample_period);

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, estimation_parameters, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(estimation_parameters, sample_period);

  /* Propagate the estimation error covariance: 3x3 matrix. */
  estimation_parameters.error_covariance =
      estimation_parameters.motion_jacobian *
          estimation_parameters.error_covariance *
          estimation_parameters.motion_jacobian.transpose() +
      estimation_parameters.process_jacobian *
          estimation_parameters.process_noise *
          estimation_parameters.process_jacobian.transpose();
}

/**
 * @brief Performs the Iterative Extended Kalman correct step.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the obejct that was
 * measured by the ego robot.
 * @param[in] robust Flag which determines whether the robust version of the
 * cost function should be used.
 */
void IEKF::correction(EstimationParameters &ego_robot,
                      const EstimationParameters &other_agent,
                      const bool robust) {

  if (robust) {
    robustCorrection(ego_robot, other_agent);
    return;
  }

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t intial_state_estimate =
      createAugmentedState(ego_robot, other_agent);

  /* Create a vector to hold the iterative state estimate. */
  augmentedState_t iterative_state_estimate = intial_state_estimate;

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_agent);

  /* Perform the iterative update.  */
  for (int i = 0; i < this->max_iterations; i++) {

    /* Calculate measurement Jacobian */
    calculateMeasurementJacobian(ego_robot, other_agent);

    /* Calculate Covariance Innovation: */
    ego_robot.innovation_covariance =
        ego_robot.measurement_jacobian * error_covariance *
            ego_robot.measurement_jacobian.transpose() +
        ego_robot.measurement_noise;

    /* Calculate Kalman Gain */
    ego_robot.kalman_gain = error_covariance *
                            ego_robot.measurement_jacobian.transpose() *
                            ego_robot.innovation_covariance.inverse();

    /* Populate the predicted measurement matrix. */
    measurement_t predicted_measurement =
        measurementModel(ego_robot, other_agent);

    /* Calculate the measurement residual. */
    ego_robot.innovation = ego_robot.measurement - predicted_measurement;

    /* Normalise the bearing residual */
    normaliseAngle(ego_robot.innovation(BEARING));

    /* Keep track of the previous estimate for the calculation of estimation
     * change at the end of the loop. */
    augmentedState_t old_estimate = iterative_state_estimate;

    ego_robot.estimation_residual =
        intial_state_estimate - iterative_state_estimate;

    normaliseAngle(ego_robot.estimation_residual(ORIENTATION));

    /* Update the iterative state estimate. */
    iterative_state_estimate =
        intial_state_estimate +
        ego_robot.kalman_gain *
            (ego_robot.innovation -
             ego_robot.measurement_jacobian * (ego_robot.estimation_residual));

    /* Break if the change between iterations converges */
    double change = (iterative_state_estimate - old_estimate).norm();
    if (change < 1e-8) {
      break;
    }
  }

  /* Resize matrices back to normal */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(ego_robot.state_estimate(ORIENTATION));

  error_covariance -= ego_robot.kalman_gain * ego_robot.innovation_covariance *
                      ego_robot.kalman_gain.transpose();

  ego_robot.error_covariance = marginalise(error_covariance);
}

/**
 * @brief A robust version of the correction function that uses the Huber cost
 * function to increase estimation error covariance of measurements that seem to
 * be outliers.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the obejct that was
 * measured by the ego robot.
 */
void IEKF::robustCorrection(EstimationParameters &ego_robot,
                            const EstimationParameters &other_agent) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t initial_state_estimate =
      createAugmentedState(ego_robot, other_agent);

  /* Create the iterative state estimate matrix. */
  augmentedState_t iterative_state_estimate = initial_state_estimate;

  /* Create and populate 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_agent);

  /* Calculate the Cholesky Decomposition of the estimation error covariance */
  Eigen::LLT<augmentedCovariance_t> error_cholesky(error_covariance);

  if (error_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the estimation error covariance");
  }

  augmentedCovariance_t error_cholesky_matrix = error_cholesky.matrixL();

  /* Calculate the Cholesky Decomposition of the sensor error covariance */
  Eigen::LLT<measurementCovariance_t> measurement_cholesky(
      ego_robot.measurement_noise);

  if (measurement_cholesky.info() != Eigen::Success) {
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the measurement error covariance");
  }

  measurementCovariance_t measurement_cholesky_matrix =
      measurement_cholesky.matrixL();

  /* Perform the iterative update.  */
  for (int i = 0; i < max_iterations; i++) {
    /* Calculate measurement Jacobian */
    calculateMeasurementJacobian(ego_robot, other_agent);

    /* Populate the predicted measurement matrix. */

    measurement_t predicted_measurement =
        measurementModel(ego_robot, other_agent);

    /* Calculate the measurement residual */
    ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

    /* Normalise the bearing residual */
    normaliseAngle(ego_robot.innovation(BEARING));

    /* Calculate the new estimation residual. */
    ego_robot.estimation_residual =
        initial_state_estimate - iterative_state_estimate;

    normaliseAngle(ego_robot.estimation_residual(ORIENTATION));

    /* Calculate the new robust estimation error covariance. */
    augmentedCovariance_t reweighted_error_covariance =
        error_cholesky_matrix *
        HuberState(ego_robot.estimation_residual, state_thresholds).inverse() *
        error_cholesky_matrix.transpose();

    /* Calculate the new robust sensor error covariance. */
    measurementCovariance_t reweighted_measurement_covariance =
        measurement_cholesky_matrix *
        HuberMeasurement(ego_robot.innovation, measurement_thresholds)
            .inverse() *
        measurement_cholesky_matrix.transpose();

    /* Calculate Covariance Innovation: */
    ego_robot.innovation_covariance =
        ego_robot.measurement_jacobian * reweighted_error_covariance *
            ego_robot.measurement_jacobian.transpose() +
        reweighted_measurement_covariance;

    /* Calculate Kalman Gain */
    ego_robot.kalman_gain = reweighted_error_covariance *
                            ego_robot.measurement_jacobian.transpose() *
                            ego_robot.innovation_covariance.inverse();

    augmentedState_t old_estimate = iterative_state_estimate;

    iterative_state_estimate =
        initial_state_estimate +
        ego_robot.kalman_gain *
            (ego_robot.innovation -
             ego_robot.measurement_jacobian * (ego_robot.estimation_residual));

    /* Break if the change between iterations converges */
    double change = (iterative_state_estimate - old_estimate).norm();
    if (change < 1e-8) {
      break;
    }
  }

  /* Resize matrices back to normal */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(ego_robot.state_estimate(ORIENTATION));

  /* Update estimation error covariance */
  error_covariance =
      (Eigen::Matrix<double, 5, 5>::Identity() -
       ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
      error_cholesky_matrix *
      HuberState(initial_state_estimate - iterative_state_estimate,
                 state_thresholds)
          .inverse() *
      error_cholesky_matrix.transpose();

  ego_robot.error_covariance = marginalise(error_covariance);
}
