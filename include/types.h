#pragma once

#include <Eigen/Dense>

namespace Filters {

/**
 * @brief The total number of system states: x coordinate [m], y-coordinate
 * [m], and orientation (heading) [rad].
 */
const unsigned short total_states{3};
/**
 * @brief The total number of system inputs: forward velocity [m/s] and
 * angular velocity [rad/s].
 */
const unsigned short total_inputs{2};
/**
 * @brief The total number of system measurements: range [m] and bearing
 * [rad].
 */
const unsigned short total_measurements{2};

/* Enums for easier accessing of variables in arrays and vectors. */
enum { FORWARD_VELOCITY = 0, ANGULAR_VELOCITY = 1 };
enum { RANGE = 0, BEARING = 1 };
enum { X = 0, Y = 1, ORIENTATION = 2 };

/* General Vector Types Definitions. */
using vector2D_t = Eigen::Matrix<double, 2, 1>;
using vector3D_t = Eigen::Matrix<double, 3, 1>;
using vector6D_t = Eigen::Matrix<double, 6, 1>;

/* General Square Matrix Type Definitions. */
using matrix2D_t = Eigen::Matrix<double, 2, 2>;
using matrix3D_t = Eigen::Matrix<double, 3, 3>;
using matrix6D_t = Eigen::Matrix<double, 6, 6>;

/* State Type Definitions. */
using state_t = vector3D_t;
using augmentedState_t = vector6D_t;

/* Covariance Type Definitions. */
using covariance_t = matrix3D_t;
using augmentedCovariance_t = matrix6D_t;

/* Information Type Definitions. */
using information_t = vector3D_t;
using augmentedInformation_t = vector6D_t;

/* Precision Type Definitions. */
using precision_t = matrix3D_t;
using augmentedPrecision_t = matrix6D_t;

/* Measurement Type Definitions. */
using measurement_t = vector2D_t;
using measurementCovariance_t = matrix2D_t;

/* Measurement Type Definitions. */
using kalmanGain_t = Eigen::Matrix<double, total_states, total_measurements>;

using augmentedKalmanGain_t =
    Eigen::Matrix<double, 2 * total_states, total_measurements>;

/* Motion Model Type Definitions. */
using motionJacobian_t = Eigen::Matrix<double, total_states, total_states>;

using processJacobian_t = Eigen::Matrix<double, total_states, total_inputs>;

using measurementJacobian_t =
    Eigen::Matrix<double, total_measurements, total_states>;

using augmentedMeasurementJacobian_t =
    Eigen::Matrix<double, total_measurements, 2 * total_states>;

using processCovariance_t = Eigen::Matrix<double, total_inputs, total_inputs>;

/* Huber Cost Function Type Definitions. */
using huberMeasurementThresholds_t =
    Eigen::Matrix<double, total_measurements, 1>;

using huberMeasurementWeights_t =
    Eigen::Matrix<double, total_measurements, total_measurements>;

using huberStateThresholds_t = Eigen::Matrix<double, 2 * total_states, 1>;
using huberStateWeights_t =
    Eigen::Matrix<double, 2 * total_states, 2 * total_states>;

} // namespace Filters
