/**
 * @file filter.h
 * @brief Header file of the parent class containing shared functionality amoung
 * cooperative localsiation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */

#ifndef INCLUDE_SRC_FILTER_H_
#define INCLUDE_SRC_FILTER_H_

#include <DataHandler.h>
#include <Eigen/Dense>

class Filter {
protected:
  /**
   * @brief The total number of system states: x coordinate [m], y-coordinate
   * [m], and orientation (heading) [rad].
   */
  static const unsigned short total_states = 3;
  /**
   * @brief The total number of system inputs: forward velocity [m/s] and
   * angular velocity [rad/s].
   */
  static const unsigned short total_inputs = 2;
  /**
   * @brief The total number of system measurements: range [m] and bearing
   * [rad].
   */
  static const unsigned short total_measurements = 2;

  /* Enums for easier accessing of variables in arrays and vectors. */
  enum { FORWARD_VELOCITY = 0, ANGULAR_VELOCITY = 1 };
  enum { RANGE = 0, BEARING = 1 };
  enum { X = 0, Y = 1, ORIENTATION = 2 };

  /* General Vector Types Definitions. */
  typedef Eigen::Matrix<double, 2, 1> vector2D_t;
  typedef Eigen::Matrix<double, total_states, 1> vector3D_t;
  typedef Eigen::Matrix<double, 2 * total_states, 1> vector6D_t;

  /* General Square Matrix Type Definitions. */
  typedef Eigen::Matrix<double, 2, 2> matrix2D_t;
  typedef Eigen::Matrix<double, total_states, total_states> matrix3D_t;
  typedef Eigen::Matrix<double, 2 * total_states, 2 * total_states> matrix6D_t;

  /* State Type Definitions. */
  typedef vector3D_t state_t;
  typedef vector6D_t augmentedState_t;

  /* Covariance Type Definitions. */
  typedef matrix3D_t covariance_t;
  typedef matrix6D_t augmentedCovariance_t;

  /* Information Type Definitions. */
  typedef vector3D_t information_t;
  typedef vector6D_t augmentedInformation_t;

  /* Precision Type Definitions. */
  typedef matrix3D_t precision_t;
  typedef matrix6D_t augmentedPrecision_t;

  /* Measurement Type Definitions. */
  typedef vector2D_t measurement_t;
  typedef matrix2D_t measurementCovariance_t;

  /* Measurement Type Definitions. */
  typedef Eigen::Matrix<double, 2 * total_states, total_measurements>
      kalmanGain_t;

  /* Motion Model Type Definitions. */
  typedef Eigen::Matrix<double, total_states, total_states> motionJacobian_t;

  typedef Eigen::Matrix<double, total_states, total_inputs> processJacobian_t;

  typedef Eigen::Matrix<double, total_measurements, 2 * total_states>
      measurementJacobian_t;

  typedef Eigen::Matrix<double, total_inputs, total_inputs> processCovariance_t;

  /* Huber Cost Function Type Definitions. */
  typedef Eigen::Matrix<double, total_measurements, 1>
      huberMeasurementThresholds_t;

  typedef Eigen::Matrix<double, total_measurements, total_measurements>
      huberMeasurementWeights_t;

  typedef Eigen::Matrix<double, 2 * total_states, 1> huberStateThresholds_t;
  typedef Eigen::Matrix<double, 2 * total_states, 2 * total_states>
      huberStateWeights_t;

  /**
   * @brief Data class housing all the data pertaining to cooperative
   * localisation (positioning).
   */
  Data::Handler &data_;

  /**
   * @struct EstimationParameters
   * @brief Houses the parameters required for performing Bayesian filtering.
   */
  struct EstimationParameters {
    unsigned short id;
    unsigned short barcode;

    /**
     * @brief Estimated robot state - x coordinate [m], y-coordinate [m], and
     * orientation (heading) [rad]: 3x1 matrix.
     * @details The state vector of the robot take the form \f[\begin{bmatrix} x
     * & y & \theta \end{bmatrix}^\top, \f] where \f$x\f$ and \f$y\f$ denotes
     * the robots 2D coordinates; and \f$\theta \f$ denotes the robots heading.
     */
    state_t state_estimate = state_t::Zero();

    /**
     * @brief Recieved range and bearing measurement: 2x1 matrix.
     */
    measurement_t measurement = measurement_t::Zero();

    /**
     * @brief The difference between the measurement recieved by the sensor and
     * the predicted measurement based on the vehicles' states.
     * @details The measurement residual is defined by the expression:
     * \f[ \tilde{\mathbf{y}}_k = \mathbf{z}_k - h(\hat{\mathbf{x}}_{k\mid
     * k-1})\f], where \f$\mathbf{z}\f$ is the measurement taken, and
     * \f$x_{k\mid k-1}\f$ is the estimated state of the robot.
     */
    measurement_t innovation = measurement_t::Zero();

    /**
     * @brief Measurement innovation covariance matrix: 2x2 matrix.
     */
    measurementCovariance_t innovation_covariance =
        measurementCovariance_t::Zero();

    /**
     * @brief The difference between the prior and posterior state estimtes.
     * @details The estimation residual is defined by the expresssion:
     * \f[\delta \mathbf{x}_k = \mathbf{x}_{k\mid k-1} - \mathbf{x}_{k\mid k}
     * \f].
     */
    augmentedState_t estimation_residual = augmentedState_t::Zero();

    /**
     * @brief Estimation Error Covariance: 3x3 matrix.
     * @details There is a high certainty in the prior value of system state,
     * therefore the prior estimation error covariance is initialised to a small
     * value.
     */
    covariance_t error_covariance = covariance_t::Identity() * 0.001;

    /**
     * @brief Kalman gain: 5x2 matrix.
     * @note The reason the Kalman gain matrix has 5 elements is because the
     * cooperative localisation (positioning) requires the estimation
     * covariances of both the ego and measured robots position [3+2]. See
     * EKF::correction.
     */
    kalmanGain_t kalman_gain = kalmanGain_t::Zero();

    /**
     * @brief Odometry process noise covariance matrix: 2x2 matrix.
     * @details The process noise covariance matrix is defined by the
     * expression:
     * \f[ w = \begin{bmatrix} q_v & 0 \\ 0 & q_\omega \end{bmatrix}, \f] where
     * \f$q_v\f$ denotes the forward velocity noise variance; and
     * \f$q_\omega\f$ denotes the angular velocity noise variance.
     * @note The process noise is assumed to be uncorrelated and therefore the
     * covariance between the forward velocity and the angular velocity is
     * assumed to be zero.
     */
    processCovariance_t process_noise = processCovariance_t::Zero();

    /**
     * @brief Measurement noise covariance matrix: 2x2 matrix.
     * @details The matrix for the measurement noise covariance matrix take the
     * form \f[ v = \begin{bmatrix} q_r & 0 \\ 0 & q_\phi \end{bmatrix}, \f]
     * where \f$q_r\f$ denotes the range noise variance; and \f$\phi_r\f$
     * denotes the bearing noise variance.
     * @note The measurement noise is assumed to be uncorrelated and therefore
     * the covariance between the range and bearing is assumed to be zero.
     */
    measurementCovariance_t measurement_noise = measurementCovariance_t::Zero();

    /**
     * @brief Jacobian matrix of the motion model evaluated in terms of the
     * systems states: x, y and orientation.
     */
    motionJacobian_t motion_jacobian = motionJacobian_t::Zero();

    /**
     * @brief Jacobian matrix of the motion model evaluated in terms of the
     * process inputs.
     */
    processJacobian_t process_jacobian = processJacobian_t::Zero();

    /**
     * @brief Jacobian of the measurement model: 2 x 6 matrix.
     */
    measurementJacobian_t measurement_jacobian = measurementJacobian_t::Zero();

    /**
     * @brief Information vector: 3x1 matrix.
     * @details Used by the information form of the (Extended) Kalman Filter.
     * The expression of the information vector take the form:
     * \f[\begin{align}\nabla &= \Sigma^{-1}\mathbf{x} \\ &= \Lambda
     * \mathbf{x}\end{align}\f], where \f$\Sigma\f$ denotes the estimation error
     * covariance matrix (Filter::EstimationParameters.error_covariance);
     * \f$\Lambda\f$ denotes the information matrix
     * (Filter::EstimationParameters.information_vector)
     * \f$\mathbf{x}\f$ denotes the state of the system
     * (Filter::EstimationParameters.state_estimate).
     */
    state_t information_vector = state_t::Zero();

    /**
     * @brief Precision matrix: 3x3 matrix.
     * @details Used by the information form of the (Extended) Kalman Filter.
     * The expression for the information matrix takes the form \f[\Lambda =
     * \Sigma^{-1}\f], where \f$\Sigma\f$ denotes the estimation error
     * covariance matrix (Filter::EstimationParameters.error_covariance).
     */
    precision_t precision_matrix = error_covariance.inverse();
  };

  /**
   * @brief Houses all estimation parameters for all robots.
   */
  std::vector<EstimationParameters> robot_parameters;

  /**
   * @brief Houses all estimation parameters for all landmarks.
   */
  std::vector<EstimationParameters> landmark_parameters;

  /**
   * @brief The thresholds for the huber measurement cost function.
   * @details The vector is of the for: range and bearing.
   */
  const huberMeasurementThresholds_t measurement_thresholds =
      (huberMeasurementThresholds_t() << 0.2, 0.01).finished();

  /**
   * @brief The thresholds for the huber state cost function.
   * @details The vector is of the for: ego x, ego y, ego orientation, agent x,
   * agent y.
   */
  const huberStateThresholds_t state_thresholds =
      (huberStateThresholds_t() << 0.15, 0.154, 0.255, 0.0104, 0.0104, 0.0)
          .finished();

  /* Filter Functionality Functions */
  virtual void prediction(const Data::Robot::Odometry &,
                          EstimationParameters &) = 0;

  virtual void correction(EstimationParameters &,
                          const EstimationParameters &) = 0;

  void motionModel(const Data::Robot::Odometry &, EstimationParameters &,
                   const double);

  void calculateMotionJacobian(const Data::Robot::Odometry &,
                               EstimationParameters &, const double);

  void calculateProcessJacobian(EstimationParameters &, const double);

  measurement_t measurementModel(EstimationParameters &,
                                 const EstimationParameters &);

  void calculateMeasurementJacobian(EstimationParameters &,
                                    const EstimationParameters &);

  matrix3D_t marginalise(const matrix6D_t &);

  state_t marginalise(const vector6D_t &, const matrix6D_t &);

  augmentedInformation_t createAugmentedVector(const state_t &,
                                               const state_t &);

  augmentedCovariance_t createAugmentedMatrix(const covariance_t &,
                                              const covariance_t &);

  void normaliseAngle(double &);

  measurement_t calculateNormalisedInnovation(const EstimationParameters &);

  augmentedState_t
  calculateNormalisedEstimationResidual(const EstimationParameters &);

  /* Huber Cost Functions. */
  huberMeasurementWeights_t
  HuberMeasurement(const measurement_t &, const huberMeasurementThresholds_t &);

  huberStateWeights_t HuberState(const augmentedState_t &,
                                 const huberStateThresholds_t &);

  matrix3D_t computePseudoInverse(const matrix3D_t &);
  matrix6D_t computePseudoInverse(const matrix6D_t &);

public:
  explicit Filter(Data::Handler &data);
  virtual ~Filter();

  void performInference();
};

#endif // INCLUDE_SRC_FILTER_H_
