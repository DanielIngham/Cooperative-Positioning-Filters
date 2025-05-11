#include "gbp.h"
#include "filter.h"

GBP::GBP(DataHandler &data, const unsigned short int window_size)
    : Filter(data), window_size_(window_size) {}

GBP::~GBP() {}

#if 0
void GBP::performInference() {
  {

    std::vector<Robot> &robots = this->data_.getRobots();
    std::vector<Landmark> landmarks = data_.getLandmarks();

    /* Loop through each timestep and perform inference.  */
    std::vector<size_t> measurement_index(data_.getNumberOfRobots(), 1);

    for (size_t k = 1; k < data_.getNumberOfSyncedDatapoints(); k++) {

      /* Perform prediction for each robot using odometry values. */
      for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

        Robot::Odometry odometry(
            robots[id].synced.odometry[k].time,
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

        for (unsigned short j = 0; j < current_measurement.subjects.size();
             j++) {
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
        }

        /* TODO: Recover the state estimate and covariance from the infomation
         * form. */
        robot_parameters[id].state_estimate =
            robot_parameters[id].precision_matrix.inverse() *
            robot_parameters[id].information_vector;

        robot_parameters[id].error_covariance =
            robot_parameters[id].precision_matrix.inverse();
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
    }
  }
}

void GBP::prediction(const Robot::Odometry &,
                     std::vector<Filter::EstimationParameters> &) {}

void GBP::correction(std::vector<EstimationParameters> &,
                     const EstimationParameters &) {}
#endif // 0
