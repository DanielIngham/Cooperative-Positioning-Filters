/**
 * @file iekf.cpp
 * @brief Implementation of the Iterative Extended Kalman Fitler implementation
 * for multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#include "iekf.h"

namespace Filters {
/**
 * @brief IEKF class constructor.
 * @details This constructor sets up the prior states and parameters to perform
 * Iterative Extended Kalman filtering.
 * @param[in] data Class containing all robot and landmark data.
 */
IEKF::IEKF(Data::Handler &data) : EKF(data) {}

/**
 * @brief Default destructor.
 */
IEKF::~IEKF() {}

/**
 * @brief Performs the Iterative Extended Kalman correct step.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the obejct that was
 * measured by the ego robot.
 */
void IEKF::correction(EstimationParameters &ego,
                      const EstimationParameters &agent) {

#ifdef ROBUST
  robustCorrection(ego_robot, other_agent);
  return;
#endif // ROBUST

  /* Create the state matrix for both robot: 5x1 matrix. */
  augmentedState_t intial_state_estimate{
      createAugmentedVector(ego.state_estimate, agent.state_estimate)};

  /* Create a vector to hold the iterative state estimate. */
  augmentedState_t iterative_state_estimate{intial_state_estimate};

  /* Create and populate new 5x5 error covariance matrix. */
  augmentedCovariance_t error_covariance{
      createAugmentedMatrix(ego.error_covariance, agent.error_covariance)};

  /* Perform the iterative update.  */
  for (int i{}; i < max_iterations_; ++i) {

    /* Calculate measurement Jacobian */
    calculateMeasurementJacobian(ego, agent);

    /* Calculate Covariance Innovation: */
    ego.innovation_covariance = ego.measurement_jacobian * error_covariance *
                                    ego.measurement_jacobian.transpose() +
                                ego.measurement_noise;

    /* Calculate Kalman Gain */
    ego.kalman_gain = error_covariance * ego.measurement_jacobian.transpose() *
                      ego.innovation_covariance.inverse();

    /* Populate the predicted measurement matrix. */
    measurement_t predicted_measurement{
        measurementModel(ego.state_estimate, agent.state_estimate)};

    /* Calculate the measurement residual. */
    ego.innovation = ego.measurement - predicted_measurement;

    /* Normalise the bearing residual */
    Data::Robot::normaliseAngle(ego.innovation(BEARING));

    /* Keep track of the previous estimate for the calculation of estimation
     * change at the end of the loop. */
    augmentedState_t old_estimate{iterative_state_estimate};

    ego.estimation_residual = intial_state_estimate - iterative_state_estimate;

    Data::Robot::normaliseAngle(ego.estimation_residual(ORIENTATION));

    /* Update the iterative state estimate. */
    iterative_state_estimate =
        intial_state_estimate +
        ego.kalman_gain * (ego.innovation - ego.measurement_jacobian *
                                                (ego.estimation_residual));

    /* Break if the change between iterations converges */
    double change{(iterative_state_estimate - old_estimate).norm()};
    static constexpr double convergence_threshold{1e-8};
    if (change < convergence_threshold) {
      break;
    }
  }

  /* Resize matrices back to normal */
  ego.state_estimate = iterative_state_estimate.head<total_states>();

  error_covariance -=
      ego.kalman_gain * ego.innovation_covariance * ego.kalman_gain.transpose();

  ego.error_covariance = error_covariance.topLeftCorner<3, 3>();
}

/**
 * @brief A robust version of the correction function that uses the Huber cost
 * function to increase estimation error covariance of measurements that seem to
 * be outliers.
 * @param[in,out] ego_robot The estimation parameters of the ego robot.
 * @param[in] other_agent The estimation parameters of the obejct that was
 * measured by the ego robot.
 */
// void IEKF::robustCorrection(EstimationParameters &ego_robot,
//                             const EstimationParameters &other_agent) {
//
//   /* Create the state matrix for both robot: 5x1 matrix. */
//   augmentedState_t initial_state_estimate{createAugmentedVector(
//       ego_robot.state_estimate, other_agent.state_estimate)};
//
//   /* Create the iterative state estimate matrix. */
//   augmentedState_t iterative_state_estimate{initial_state_estimate};
//
//   /* Create and populate 5x5 error covariance matrix. */
//   augmentedCovariance_t error_covariance{createAugmentedMatrix(
//       ego_robot.error_covariance, other_agent.error_covariance)};
//
//   /* Calculate the Cholesky Decomposition of the estimation error covariance
//   */ Eigen::LLT<augmentedCovariance_t> error_cholesky{error_covariance};
//
//   if (error_cholesky.info() != Eigen::Success) {
//     throw std::runtime_error(
//         "An error has occurred with calculating the Cholesky decomposition of
//         " "the estimation error covariance");
//   }
//
//   augmentedCovariance_t error_cholesky_matrix{error_cholesky.matrixL()};
//
//   /* Calculate the Cholesky Decomposition of the sensor error covariance */
//   Eigen::LLT<measurementCovariance_t> measurement_cholesky(
//       ego_robot.measurement_noise);
//
//   if (measurement_cholesky.info() != Eigen::Success) {
//     throw std::runtime_error(
//         "An error has occurred with calculating the Cholesky decomposition of
//         " "the measurement error covariance");
//   }
//
//   measurementCovariance_t measurement_cholesky_matrix =
//       measurement_cholesky.matrixL();
//
//   /* Perform the iterative update.  */
//   for (int i{}; i < max_iterations_; i++) {
//     /* Calculate measurement Jacobian */
//     calculateMeasurementJacobian(ego_robot, other_agent);
//
//     /* Populate the predicted measurement matrix. */
//
//     measurement_t predicted_measurement{
//         measurementModel(ego_robot, other_agent)};
//
//     /* Calculate the measurement residual */
//     ego_robot.innovation = (ego_robot.measurement - predicted_measurement);
//
//     /* Normalise the bearing residual */
//     Data::Robot::normaliseAngle(ego_robot.innovation(BEARING));
//
//     /* Calculate the new estimation residual. */
//     ego_robot.estimation_residual =
//         initial_state_estimate - iterative_state_estimate;
//
//     Data::Robot::normaliseAngle(ego_robot.estimation_residual(ORIENTATION));
//
//     /* Calculate the new robust estimation error covariance. */
//     augmentedCovariance_t reweighted_error_covariance{
//         error_cholesky_matrix *
//         HuberState(ego_robot.estimation_residual, state_thresholds).inverse()
//         * error_cholesky_matrix.transpose()};
//
//     /* Calculate the new robust sensor error covariance. */
//     measurementCovariance_t reweighted_measurement_covariance{
//         measurement_cholesky_matrix *
//         HuberMeasurement(ego_robot.innovation, measurement_thresholds)
//             .inverse() *
//         measurement_cholesky_matrix.transpose()};
//
//     /* Calculate Covariance Innovation: */
//     ego_robot.innovation_covariance =
//         ego_robot.measurement_jacobian * reweighted_error_covariance *
//             ego_robot.measurement_jacobian.transpose() +
//         reweighted_measurement_covariance;
//
//     /* Calculate Kalman Gain */
//     ego_robot.kalman_gain = reweighted_error_covariance *
//                             ego_robot.measurement_jacobian.transpose() *
//                             ego_robot.innovation_covariance.inverse();
//
//     augmentedState_t old_estimate{iterative_state_estimate};
//
//     iterative_state_estimate =
//         initial_state_estimate +
//         ego_robot.kalman_gain *
//             (ego_robot.innovation -
//              ego_robot.measurement_jacobian *
//              (ego_robot.estimation_residual));
//
//     /* Break if the change between iterations converges */
//     double change{(iterative_state_estimate - old_estimate).norm()};
//
//     static constexpr double convergence_threshold{1e-8};
//     if (change < convergence_threshold) {
//       break;
//     }
//   }
//
//   /* Resize matrices back to normal */
//   ego_robot.state_estimate = iterative_state_estimate.head<total_states>();
//
//   /* Update estimation error covariance */
//   error_covariance =
//       (augmentedCovariance_t::Identity() -
//        ego_robot.kalman_gain * ego_robot.measurement_jacobian) *
//       error_cholesky_matrix *
//       HuberState(initial_state_estimate - iterative_state_estimate,
//                  state_thresholds)
//           .inverse() *
//       error_cholesky_matrix.transpose();
//
//   ego_robot.error_covariance = error_covariance.topLeftCorner<3, 3>();
// }

} // namespace Filters
