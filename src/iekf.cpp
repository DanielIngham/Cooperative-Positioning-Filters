/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "iekf.h"
#include "filter.h"
#include <DataHandler/Landmark.h>
#include <DataHandler/Robot.h>
#include <Eigen/Cholesky>
#include <cmath>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>

/**
 * @brief EKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Extended Kalman filtering.
 * @param[in] data Class containing all robot data.
 */
IEKF::IEKF(DataHandler &data) : Filter(data) {}

/**
 * @brief Default destructor.
 */
IEKF::~IEKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void IEKF::performInference() {
  std::vector<Robot> &robots = this->data_.getRobots();
  std::vector<Landmark> landmarks = data_.getLandmarks();

  /* Loop through each timestep and perform inference.  */
  std::vector<size_t> measurement_index(data_.getNumberOfRobots(), 0);

  std::ofstream file("output/normalised_innovation.dat");
  if (!file.is_open()) {
    throw std::runtime_error("Unable to create normalised_innovation.dat");
  }
  std::ofstream file2("output/normalised_estimation.dat");

  if (!file2.is_open()) {
    throw std::runtime_error("Unable to create normalised_innovation.dat");
  }

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

      /* Check if a measurement is available. */
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

        bool robust = true;

        if (robust) {

          Eigen::Matrix<double, total_measurements, 1> measurement_tau;
          measurement_tau << 0.2, 0.01;

          Eigen::Matrix<double, 2 + total_states, 1> state_tau;
          state_tau << 0.155, 0.154, 0.255, 10.0051104, 10.001104;

          robustCorrection(robot_parameters[id], measured_object,
                           measurement_tau, state_tau);

        } else {

          correction(robot_parameters[id], measured_object);
        }

        Eigen::Matrix<double, total_measurements, 1> measurement =
            calculateNormalisedMeasurementResidual(robot_parameters[id]);

        file << measurement(RANGE) << '\t' << measurement(BEARING) << '\n';

        Eigen::Matrix<double, 2 + total_states, 1> estimation =
            calculateNormalisedEstimationResidual(robot_parameters[id]);

        file2 << estimation(X) << '\t' << estimation(Y) << '\t'
              << estimation(ORIENTATION) << '\t' << estimation(X + 2) << '\t'
              << estimation(Y + 2) << '\n';

        /* Update the robot state data structure. */
        robots[id].synced.states[k].x = robot_parameters[id].state_estimate(X);
        robots[id].synced.states[k].y = robot_parameters[id].state_estimate(Y);
        robots[id].synced.states[k].orientation =
            robot_parameters[id].state_estimate(ORIENTATION);

        double error = robots[id].synced.states[k].orientation -
                       robots[id].groundtruth.states[k].orientation;

        normaliseAngle(error);
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
 *
 * @details The motion model used for the extended kalman filter prediction take
 * the form
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
 * EKF::EstimationParameters::process_noise).
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
 * @brief Performs the Extended Kalman correct step.
 * @param[in,out] ego_robot The parameters required by the Extended
 * Kalman filter to perform the correction step.
 * @param[in] other_object The robot that was measured by the ego robot.
 */
void IEKF::correction(EstimationParameters &ego_robot,
                      const EstimationParameters &other_object) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t intial_state_estimate =
      createAugmentedState(ego_robot, other_object);

  /* Create a vector to hold the iterative state estimate. */
  augmentedState_t iterative_state_estimate = intial_state_estimate;

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_object);

  /* Perform the iterative update.  */
  for (int i = 0; i < this->max_iterations; i++) {

    /* Calculate measurement Jacobian */
    calculateMeasurementJacobian(ego_robot, other_object);

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
        measurementModel(ego_robot, other_object);

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

void IEKF::robustCorrection(EstimationParameters &ego_robot,
                            const EstimationParameters &other_object,
                            huberMeasurementThresholds_t &measurement_taus,
                            huberStateThresholds_t &state_taus) {

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t initial_state_estimate =
      createAugmentedState(ego_robot, other_object);

  /* Create the iterative state estimate matrix. */
  augmentedState_t iterative_state_estimate = initial_state_estimate;

  /* Create and populate 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance =
      createAugmentedCovariance(ego_robot, other_object);

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

    /* Calculate the new estimation residual. */
    ego_robot.estimation_residual =
        initial_state_estimate - iterative_state_estimate;

    normaliseAngle(ego_robot.estimation_residual(ORIENTATION));

    /* Calculate the new robust estimation error covariance. */
    augmentedCovariance_t reweighted_error_covariance =
        error_cholesky_matrix *
        HuberState(ego_robot.estimation_residual, state_taus).inverse() *
        error_cholesky_matrix.transpose();

    /* Calculate the new robust sensor error covariance. */
    measurementCovariance_t reweighted_measurement_covariance =
        measurement_cholesky_matrix *
        HuberMeasurement(ego_robot.innovation, measurement_taus).inverse() *
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
  /* NOTE: The resize means all information regarding the observed robots
   * states is lost. This operates on the assumptoin that the error covariance
   * is uncorrelated, which may lead to the estimator being over confident in
   * bad estimates. */
  ego_robot.state_estimate = iterative_state_estimate.head<total_states>();

  /* Normalise the orientation estimate between -180 and 180. */
  normaliseAngle(ego_robot.state_estimate(ORIENTATION));

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix by incorporating the
   * marginalising the contributions of the error covariance from the other
   * robot states into the covariance of the ego robot. */
  /* Update estimation error covariance */
  error_covariance =
      (Eigen::Matrix<double, 5, 5>::Identity() -
       ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
      error_cholesky_matrix *
      HuberState(initial_state_estimate - iterative_state_estimate, state_taus)
          .inverse() *
      error_cholesky_matrix.transpose();

  ego_robot.error_covariance = marginalise(error_covariance);
}
