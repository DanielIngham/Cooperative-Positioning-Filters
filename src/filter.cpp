#include "filter.h"
#include <DataHandler/DataHandler.h>

/**
 * @brief Assigns fields data based on datahandler input.
 * @param[in] data Class containing all robot data.
 */
Filter::Filter(DataHandler &data) : data_(data) {

  std::vector<Robot> &robots = data_.getRobots();

  /* Populate the Estimation parameters for each robot. */
  for (unsigned short id = 0; id < data_.getNumberOfRobots(); id++) {

    EstimationParameters initial_parameters;

    /* Assume known prior. This is done by setting the first value of the
     * estimated values to the groundtruth. */
    robots[id].synced.states.push_back(
        Robot::State(robots[id].groundtruth.states.front()));

    /* Initial state: 3x1 Matrix. */
    initial_parameters.state_estimate << robots[id].synced.states.front().x,
        robots[id].synced.states.front().y,
        robots[id].synced.states.front().orientation;

    /* Populate odometry error covariance matrix: 2x2 matrix. */
    initial_parameters.process_noise.diagonal().topRows(total_inputs)
        << robots[id].forward_velocity_error.variance,
        robots[id].angular_velocity_error.variance;

    /* Populate measurement error covariance matrix: 2x2 matrix. */
    initial_parameters.measurement_noise.diagonal().topRows(total_measurements)
        << robots[id].range_error.variance,
        robots[id].bearing_error.variance;

    robot_parameters.push_back(initial_parameters);
  }

  /* Populate the estimation parameters for each landmark. */
  std::vector<Landmark> landmarks = data_.getLandmarks();

  for (unsigned short id = 0; id < data_.getNumberOfLandmarks(); id++) {
    EstimationParameters initial_parameters;

    initial_parameters.state_estimate << landmarks[id].x, landmarks[id].y, 0.0;

    /* The landmark only has two states: x and y coordintate.
     * NOTE: Although the landmark only has two states, the same data structure
     * is used for both the robot and landmark for compatibility with the
     * EFK::correction function, since only the x and y coordinate and thier
     * corresponding error covariances are used for the measurement update. */
    initial_parameters.error_covariance.diagonal().topRows(total_states - 1)
        << landmarks[id].x_std_dev * landmarks[id].x_std_dev,
        landmarks[id].y_std_dev * landmarks[id].y_std_dev;

    landmark_parameters.push_back(initial_parameters);
  }
}

Filter::~Filter() = default;
