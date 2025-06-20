/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "ekf.h"
#include "filter.h"
#include <DataHandler/Landmark.h>
#include <DataHandler/Robot.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

/**
 * @brief EKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Kalman filtering.
 * @param[in] data Class containing all robot data.
 */
EKF::EKF(DataHandler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
EKF::~EKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void EKF::performInference() {
  std::vector<Robot> &robots = this->data_.getRobots();
  std::vector<Landmark> landmarks = data_.getLandmarks();

  /* Loop through each timestep and perform inference.  */
  std::vector<size_t> measurement_index(data_.getNumberOfRobots(), 0);

  for (size_t k = 1; k < data_.getNumberOfSyncedDatapoints(); k++) {

    /* Perform prediction for each robot using odometry values. */
    for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

      Robot::Odometry odometry(robots[id].synced.odometry[k].time,
                               robots[id].synced.odometry[k].forward_velocity,
                               robots[id].synced.odometry[k].angular_velocity);

      prediction(odometry, robot_parameters[id]);

      /* Update the robot state data structure. */
      robots[id].synced.states.push_back(
          Robot::State(robots[id].groundtruth.states[k].time,
                       robot_parameters[id].state_estimate(X),
                       robot_parameters[id].state_estimate(Y),
                       robot_parameters[id].state_estimate(ORIENTATION)));
    }

    /* If a measurements are available, loop through each measurement
     * and update the estimate. */
    for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

      /* Range check. */
      if (measurement_index[id] >= robots[id].synced.measurements.size()) {
        continue;
      }

      /* Check if a measurement is not available for this time stamp, then skip.
       */
      if (std::round(
              (robots[id].synced.measurements[measurement_index[id]].time -
               robots[id].synced.odometry[k].time) *
              10000.0) /
              10000.0 !=
          0.0) {
        continue;
      }

      /* Loop through the measurements taken and perform the measurement
       * update for each robot.
       * NOTE: This operation uses the assumption that the measurements fo the
       * indpendent robots/landmarks are independent of one another.
       */
      const Robot::Measurement &current_measurement =
          robots[id].synced.measurements[measurement_index[id]];

      for (unsigned short j = 0; j < current_measurement.subjects.size(); j++) {
        /* Find the subject for whom the barcode belongs to. */
        int subject_id = data_.getID(current_measurement.subjects[j]);

        if (-1 == subject_id) {
          continue;
        }
        /* Populate the measurement matrix required for the correction step.
         * Remove any noise bias from the measurement.
         */
        robot_parameters[id].measurement << current_measurement.ranges[j],
            current_measurement.bearings[j];

        /* The datahandler first assigns the ID to the robots then the
         * landmarks. Therefore if the ID is less than or equal to the number
         * of robots, then it belongs to a robot, otherwise it belong to a
         * landmark. */
        EstimationParameters measured_object;

        if (subject_id <= data_.getNumberOfRobots()) {
          unsigned short index = subject_id - 1;
          measured_object = robot_parameters[index];

        } else {
          unsigned short index = subject_id - data_.getNumberOfRobots() - 1;
          measured_object = landmark_parameters[index];
        }

        total_observations++;

        /* Determine whether the filter should use the huber cost function for
         * the correction. */
        bool robust = true;
        if (robust) {
          huberMeasurementThresholds_t measurement_tau =
              (huberMeasurementThresholds_t() << 0.2, 0.01).finished();

          huberStateThresholds_t state_tau =
              (huberStateThresholds_t() << 0.15, 0.154, 0.255, 0.0104, 0.0104)
                  .finished();

          robustCorrection(robot_parameters[id], measured_object,
                           measurement_tau, state_tau);
        } else {
          correction(robot_parameters[id], measured_object);
        }

        /* Update the robot state data structure. */
        robots[id].synced.states[k].x = robot_parameters[id].state_estimate(X);
        robots[id].synced.states[k].y = robot_parameters[id].state_estimate(Y);
        robots[id].synced.states[k].orientation =
            robot_parameters[id].state_estimate(ORIENTATION);
      }
      measurement_index[id] += 1;
    }
  }

  /* Calculate the inference error. */
  for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {
    robots[id].calculateStateError();
  }
}

/**
 * @brief performs the prediction step of the Extended Kalman filter.
 * @param[in] odometry The prior inputs into the system comprising a forward and
 * angular velocity.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the prediction step.
 */
void EKF::prediction(const Robot::Odometry &odometry,
                     EstimationParameters &estimation_parameters) {

  double sample_period = data_.getSamplePeriod();

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
 * @brief Performs the Extended Kalman correct step.
 * @param[in,out] estimation_parameters The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_robot The robot that was measured by the ego robot.
 * @param[in] update Determines whether the state and covariance should be
 * updated.
 * @details The measusurement model for the measurement taken from ego vehicle
 * \f$i\f$ to vehicle \f$j\f$ used for the correction step takes the form
 * \f[ \begin{bmatrix} r_{ij}^{(t)} \\ \phi_{ij}^{(t)}\end{bmatrix} =
 * \begin{bmatrix}\sqrt{(x_j^{(t)} - x_i^{(t)})^2 + (y_j^{(t)} - y_i^{(t)})^2} +
 * q_r \\ \text{atan2}\left(\frac{y_j^{(t)}-y_i^{(t)}}{x_j^{(t)}-x_i^{(t)}
 * }\right) - \theta_i^{(t)} + q_\phi\end{bmatrix}, \f] where \f$x\f$ and
 * \f$y\f$ denote the robots coordinates; \f$\theta\f$ denotes the ego robots
 * orientation (heading); and \f$q_r\f$ and \f$q_\omega\f$ denote the Gaussian
 * distributed measurement noise (See
 * Filter::EstimationParameters.measurement_noise).
 *
 * @note Cooperative Localisation (Positioning) involves robots that share thier
 * state and estimation error covariances when one robot measures the other. As
 * a result, the estimation error covariance needs to be augmented from a 3x3 to
 * a 5x5 matrix to house the error covariance of both the ego vehicle (\f$i\f$)
 * and the measured vehicle (\f$j\f$):
 * \f[\mathbf{P} = \begin{bmatrix} \mathbf{P}_i & \mathbf{0} \\ \mathbf{0} &
 * \mathbf{P}_j \end{bmatrix}, \f] where \f$\mathbf{P}_i\f$ and
 * \f$\mathbf{P}_j\f$ are the estimation error covariance of the ego robot
 * \f$i\f$ and the observed robot \f$j\f$ respectively.
 */
void EKF::correction(EstimationParameters &ego_robot,
                     const EstimationParameters &other_object) {

  /* Calculate measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_object);

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_object);

  /* Calculate Covariance Innovation: */
  ego_robot.innovation_covariance =
      ego_robot.measurement_jacobian * error_covariance *
          ego_robot.measurement_jacobian.transpose() +
      ego_robot.measurement_noise;

  /* Calculate Kalman Gain */
  ego_robot.kalman_gain = error_covariance *
                          ego_robot.measurement_jacobian.transpose() *
                          ego_robot.innovation_covariance.inverse();

  /* Update estimation error covariance */
  error_covariance -= ego_robot.kalman_gain * ego_robot.innovation_covariance *
                      ego_robot.kalman_gain.transpose();

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t state_estimate =
      createAugmentedState(ego_robot, other_object);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_object);

  /* Calculate the innovation: the difference between the measurement
   * and the predicted measurement based on the estimated states of both robots.
   */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual. */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Update the state using the measurement. */
  state_estimate += ego_robot.kalman_gain * ego_robot.innovation;

  normaliseAngle(state_estimate(ORIENTATION));

  /* Resize matrices back to normal */
  ego_robot.state_estimate = state_estimate.head<total_states>();

  ego_robot.error_covariance = marginalise(error_covariance);
}

void EKF::robustCorrection(EstimationParameters &ego_robot,
                           const EstimationParameters &other_object,
                           huberMeasurementThresholds_t &measurement_tau,
                           huberStateThresholds_t &state_tau) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t intial_state_estimate =
      createAugmentedState(ego_robot, other_object);

  augmentedState_t iterative_state_estimate = intial_state_estimate;

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_object);

  /* Calculate the Cholesky Decomposition of the estimation error covariance */
  Eigen::LLT<augmentedCovariance_t> error_cholesky(error_covariance);

  if (error_cholesky.info() != Eigen::Success) {
    throw std::runtime_error("[1] An error has occurred with calculating the "
                             "Cholesky decomposition of "
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

  /* Calculate measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_object);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_object);

  /* Calculate the measurement residual: the difference between the
   * measurement and the calculate measurement based on the estimated states
   * of both robots. */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the bearing residual */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Calculate the new robust sensor error covariance. */
  measurementCovariance_t reweighted_measurement_covariance =
      measurement_cholesky_matrix *
      HuberMeasurement(ego_robot.innovation, measurement_tau).inverse() *
      measurement_cholesky_matrix.transpose();

  /* Calculate Covariance Innovation. */
  ego_robot.innovation_covariance =
      ego_robot.measurement_jacobian * error_covariance *
          ego_robot.measurement_jacobian.transpose() +
      reweighted_measurement_covariance;

  /* Calculate Kalman Gain. */
  ego_robot.kalman_gain = error_covariance *
                          ego_robot.measurement_jacobian.transpose() *
                          ego_robot.innovation_covariance.inverse();

  /* Update the inititial state estimate. */
  iterative_state_estimate =
      intial_state_estimate + ego_robot.kalman_gain * ego_robot.innovation;

  normaliseAngle(iterative_state_estimate(ORIENTATION));

  ego_robot.estimation_residual =
      intial_state_estimate - iterative_state_estimate;

  normaliseAngle(ego_robot.estimation_residual(ORIENTATION));

  /* Resize matrices back to normal */
  /* NOTE: The resize means all information regarding the observed robots
   * states is lost. This operates on the assumption that the error covariance
   * is uncorrelated, which may lead to the estimator being over confident in
   * bad estimates. */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Calculate the reweighted error covariance. */
  error_covariance =
      (augmentedCovariance_t::Identity() -
       ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
      error_cholesky_matrix *
      HuberState(ego_robot.estimation_residual, state_tau).inverse() *
      error_cholesky_matrix.transpose();

  /* Marginalise the 5x5 back to a 3x3. */
  ego_robot.error_covariance = marginalise(error_covariance);
}
