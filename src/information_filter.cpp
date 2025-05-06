#include "information_filter.h"

InformationFilter::InformationFilter(DataHandler &data) : Filter(data) {}

InformationFilter::~InformationFilter() {}

void InformationFilter::performInference() {

  std::vector<Robot> &robots = this->data_.getRobots();
  std::vector<Landmark> landmarks = data_.getLandmarks();

  /* Loop through each timestep and perform inference.  */
  std::vector<size_t> measurement_index(data_.getNumberOfRobots(), 1);

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

#if 1
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

        correction(robot_parameters[id], measured_object);
        robot_parameters[id].state_estimate =
            robot_parameters[id].information_matrix.inverse() *
            robot_parameters[id].information_vector;
      }
      /* TODO: Recover the state estimate and covariance from the infomation
       * form. */

      robot_parameters[id].error_covariance =
          robot_parameters[id].information_matrix.inverse();
      /* Normalise the orientation estimate between -180 and 180. */

      while (robot_parameters[id].state_estimate(ORIENTATION) >= M_PI)
        robot_parameters[id].state_estimate(ORIENTATION) -= 2.0 * M_PI;

      while (robot_parameters[id].state_estimate(ORIENTATION) < -M_PI)
        robot_parameters[id].state_estimate(ORIENTATION) += 2.0 * M_PI;

      /* Update the robot state data structure. */
      robots[id].synced.states[k].x = robot_parameters[id].state_estimate(X);
      robots[id].synced.states[k].y = robot_parameters[id].state_estimate(Y);
      robots[id].synced.states[k].orientation =
          robot_parameters[id].state_estimate(ORIENTATION);
      measurement_index[id] += 1;
    }
#endif // 0
  }
}

void InformationFilter::prediction(const Robot::Odometry &odometry,
                                   EstimationParameters &ego_robot) {

  double sample_period = data_.getSamplePeriod();

  /* Make the prediction using the motion model: 3x1 matrix. */
  ego_robot.state_estimate << ego_robot.state_estimate(X) +
                                  odometry.forward_velocity * sample_period *
                                      std::cos(ego_robot.state_estimate(
                                          ORIENTATION)),
      ego_robot.state_estimate(Y) +
          odometry.forward_velocity * sample_period *
              std::sin(ego_robot.state_estimate(ORIENTATION)),
      ego_robot.state_estimate(ORIENTATION) +
          odometry.angular_velocity * sample_period;

  /* Normalise the orientation estimate between -180 and 180. */
  while (ego_robot.state_estimate(ORIENTATION) >= M_PI)
    ego_robot.state_estimate(ORIENTATION) -= 2.0 * M_PI;

  while (ego_robot.state_estimate(ORIENTATION) < -M_PI)
    ego_robot.state_estimate(ORIENTATION) += 2.0 * M_PI;

  /* Calculate the Motion Jacobian: 3x3 matrix. */
  ego_robot.motion_jacobian << 1, 0,
      -odometry.forward_velocity * sample_period *
          std::sin(ego_robot.state_estimate(ORIENTATION)),
      0, 1,
      odometry.forward_velocity * sample_period *
          std::cos(ego_robot.state_estimate(ORIENTATION)),
      0, 0, 1;

  /* Calculate the process noise Jacobian: 3x2 matrix. */
  ego_robot.process_jacobian
      << sample_period * std::cos(ego_robot.state_estimate(ORIENTATION)),
      0, sample_period * std::sin(ego_robot.state_estimate(ORIENTATION)), 0, 0,
      sample_period;

  /* Propagate the estimation error covariance: 3x3 matrix. */
  ego_robot.error_covariance =
      ego_robot.motion_jacobian * ego_robot.error_covariance *
          ego_robot.motion_jacobian.transpose() +
      ego_robot.process_jacobian * ego_robot.process_noise *
          ego_robot.process_jacobian.transpose();

  /* Update the information form:
   * NOTE: this is the only place where the information filter differs from the
   * Kalman filter in the prediction step.*/
  ego_robot.information_matrix = ego_robot.error_covariance.inverse();

  ego_robot.information_vector =
      ego_robot.information_matrix * ego_robot.state_estimate;
}

void InformationFilter::correction(EstimationParameters &ego_robot,
                                   const EstimationParameters &other_object) {

  /* Calculate measurement Jacobian */
  double x_difference =
      other_object.state_estimate(X) - ego_robot.state_estimate(X);

  double y_difference =
      other_object.state_estimate(Y) - ego_robot.state_estimate(Y);

  double denominator =
      std::sqrt(x_difference * x_difference + y_difference * y_difference);

  const double MIN_DISTANCE = 1e-6;
  if (denominator < MIN_DISTANCE) {
    denominator = MIN_DISTANCE;
  }

  Eigen::Matrix<double, total_measurements, 2 * total_states>
      measurement_jacobian;

  measurement_jacobian << -x_difference / denominator,
      -y_difference / denominator, 0, x_difference / denominator,
      y_difference / denominator, 0, y_difference / (denominator * denominator),
      -x_difference / (denominator * denominator), -1,
      -y_difference / (denominator * denominator),
      x_difference / (denominator * denominator), 0;
  /* NOTE: Measurement noise Jacobian is identity. No need to calculate. */

  /* Populate the predicted measurement matrix. */
  Eigen::Matrix<double, total_measurements, 1> predicted_measurement;

  predicted_measurement << std::sqrt((x_difference * x_difference) +
                                     (y_difference * y_difference)),
      std::atan2(y_difference, x_difference) -
          ego_robot.state_estimate[ORIENTATION];

  /* Calculate the measurement residual: the difference between the measurement
   * and the calculate measurement based on the estimated states of both robots.
   */
  Eigen::Matrix<double, total_measurements, 1> measurement_residual =
      (ego_robot.measurement - predicted_measurement);

  /* Normalise the angle residual */
  while (measurement_residual(BEARING) >= M_PI)
    measurement_residual(BEARING) -= 2.0 * M_PI;

  while (measurement_residual(BEARING) < -M_PI)
    measurement_residual(BEARING) += 2.0 * M_PI;

  /* Create the state matrix for both robot: 5x1 matrix. */
  Eigen::Matrix<double, 2 * total_states, 1> estimated_state;
  estimated_state.setZero();

  estimated_state.head<total_states>() = ego_robot.state_estimate;
  estimated_state.tail<total_states>() =
      other_object.state_estimate.head<total_states>();

  /* Calculate the Information contribution */
  Eigen::Matrix<double, 2 * total_states, 2 * total_states>
      information_matrix_contribution =
          measurement_jacobian.transpose() *
          ego_robot.measurement_noise.inverse() * measurement_jacobian;

  Eigen::Matrix<double, 2 * total_states, 1> information_vector_contribution =
      measurement_jacobian.transpose() * ego_robot.measurement_noise.inverse() *
      (measurement_residual + measurement_jacobian * estimated_state);

  /* Create a temporary augmented matrix  containing the information matrix of
   * both objects. */
  Eigen::Matrix<double, 2 * total_states, 2 * total_states> information_matrix;

  information_matrix.topLeftCorner<total_states, total_states>() =
      ego_robot.information_matrix;

  information_matrix.bottomRightCorner<total_states, total_states>() =
      other_object.information_matrix
          .topLeftCorner<total_states, total_states>();

  /* Create a temporary augmented vector containing the information vector of
   * both objects. */
  Eigen::Matrix<double, 2 * total_states, 1> information_vector;

  information_vector.head<total_states>() =
      ego_robot.information_vector.head<total_states>();

  information_vector.tail<total_states>() =
      other_object.information_vector.head<total_states>();

  /* Add the information contribution. */
  information_matrix += information_matrix_contribution;
  information_vector += information_vector_contribution;

  /* Schur complement-based error covariance marginalisation. This is used to
   * marginalise the 5x5 matrix to a 3x3 matrix */
  ego_robot.information_matrix =
      information_matrix.topLeftCorner<total_states, total_states>() -
      information_matrix.topRightCorner<total_states, total_states>() *
          information_matrix.bottomRightCorner<total_states, total_states>()
              .inverse() *
          information_matrix.bottomLeftCorner<total_states, total_states>();

  ego_robot.information_vector =
      information_vector.head<total_states>() -
      information_matrix.topRightCorner<total_states, total_states>() *
          information_matrix.bottomRightCorner<total_states, total_states>()
              .inverse() *
          information_vector.tail<total_states>();
}
