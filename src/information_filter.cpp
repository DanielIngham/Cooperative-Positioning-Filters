#include "information_filter.h"

InformationFilter::InformationFilter(DataHandler &data) : Filter(data) {}

InformationFilter::~InformationFilter() {}

void InformationFilter::performInference() {

  std::vector<Robot> &robots = this->data_.getRobots();

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

    /* If measurements are available, loop through each measurement
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

        correction(robot_parameters[id], measured_object);

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

void InformationFilter::prediction(const Robot::Odometry &odometry,
                                   EstimationParameters &ego_robot) {

  const double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model: 3x1 matrix. */
  motionModel(odometry, ego_robot, sample_period);

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  calculateMotionJacobian(odometry, ego_robot, sample_period);

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  calculateProcessJacobian(ego_robot, sample_period);

  /* Propagate the estimation information: 3x3 matrix. */
  ego_robot.precision_matrix =
      (ego_robot.motion_jacobian * ego_robot.precision_matrix.inverse() *
           ego_robot.motion_jacobian.transpose() +
       ego_robot.process_jacobian * ego_robot.process_noise *
           ego_robot.process_jacobian.transpose())
          .inverse();

  ego_robot.information_vector =
      ego_robot.precision_matrix * ego_robot.state_estimate;
}

void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_object) {

  /* Calculate measurement Jacobian */
  calculateMeasurementJacobian(ego_robot, other_object);

  /* Populate the predicted measurement matrix. */
  measurement_t predicted_measurement =
      measurementModel(ego_robot, other_object);

  /* Calculate the measurement residual: the difference between the measurement
   * and the calculate measurement based on the estimated states of both robots.
   */
  ego_robot.innovation = (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual */
  normaliseAngle(ego_robot.innovation(BEARING));

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t estimated_state =
      createAugmentedState(ego_robot, other_object);

  /* Calculate the Information contribution */
  augmentedPrecision_t precision_matrix_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() * ego_robot.measurement_jacobian;

  augmentedState_t information_vector_contribution =
      ego_robot.measurement_jacobian.transpose() *
      ego_robot.measurement_noise.inverse() *
      (ego_robot.innovation + ego_robot.measurement_jacobian * estimated_state);

  /* Create a temporary augmented matrix  containing the information matrix of
   * both objects. */
  augmentedPrecision_t precision_matrix =
      createAugmentedPrecision(ego_robot, other_object);

  /* Create a temporary augmented vector containing the information vector of
   * both objects. */
  augmentedState_t information_vector = precision_matrix * estimated_state;

  /* Add the information contribution. */
  precision_matrix += precision_matrix_contribution;
  information_vector += information_vector_contribution;

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix */
  ego_robot.precision_matrix = marginalise(precision_matrix);

  ego_robot.information_vector =
      marginalise(information_vector, precision_matrix);

  ego_robot.state_estimate =
      ego_robot.precision_matrix.inverse() * ego_robot.information_vector;
}
