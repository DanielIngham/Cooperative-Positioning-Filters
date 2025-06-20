#ifndef INCLUDE_SRC_FILTER_H_
#define INCLUDE_SRC_FILTER_H_

#include <DataHandler/DataHandler.h>
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

  enum { FORWARD_VELOCITY = 0, ANGULAR_VELOCITY = 1 };
  enum { RANGE = 0, BEARING = 1 };
  enum { X = 0, Y = 1, ORIENTATION = 2 };

  typedef Eigen::Matrix<double, total_states, 1> state_t;
  typedef Eigen::Matrix<double, 2 + total_states, 1> augmentedState_t;

  typedef Eigen::Matrix<double, total_states, total_states> matrix3D_t;
  typedef Eigen::Matrix<double, 2 + total_states, 2 + total_states> matrix5D_t;

  typedef matrix3D_t covariance_t;
  typedef matrix5D_t augmentedCovariance_t;

  typedef matrix5D_t augmentedPrecision_t;

  typedef Eigen::Matrix<double, total_measurements, 1> measurement_t;
  typedef Eigen::Matrix<double, total_measurements, total_measurements>
      measurementCovariance_t;

  typedef Eigen::Matrix<double, 2 + total_states, total_measurements>
      kalmanGain_t;

  typedef Eigen::Matrix<double, total_states, total_states> motionJacobian_t;
  typedef Eigen::Matrix<double, total_states, total_inputs> processJacobian_t;
  typedef Eigen::Matrix<double, total_measurements, 2 + total_states>
      measurementJacobian_t;

  typedef Eigen::Matrix<double, total_inputs, total_inputs> processCovariance_t;

  typedef Eigen::Matrix<double, total_measurements, 1>
      huberMeasurementThresholds_t;
  typedef Eigen::Matrix<double, total_measurements, total_measurements>
      huberMeasurementWeights_t;

  typedef Eigen::Matrix<double, 2 + total_states, 1> huberStateThresholds_t;
  typedef Eigen::Matrix<double, 2 + total_states, 2 + total_states>
      huberStateWeights_t;

  /**
   * @brief Data class housing all the data pertaining to cooperative
   * localisation (positioning).
   */
  DataHandler &data_;

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
     * @brief Jacobian matrix of the motion model.
     * @details The formula used for the calculation of the motion model
     * Jacobian takes the form:
     * \f[ F = \begin{bmatrix} 1 & 0 & -\tilde{v}\Delta t \sin(\theta) \\ 0 & 1
     * & \tilde{v} \Delta t \cos(\theta) \\ 0 & 0 & 1 \end{bmatrix}, \f]
     * where \f$\theta\f$ denotes the heading (orientation) of the ego vehicle;
     * and \f$\tilde{v}\f$ denotes the forward velocity. The forward velocity is
     * a random variable with Gaussian distributed noise \f$\mathcal{N}(0,w)\f$
     * , where \f$w\f$ is defined by the covariance matrix
     * EKF::EstimationParameters.measurement_noise. See EKF::prediction for
     * information on the motion model from which this was derived.
     */
    motionJacobian_t motion_jacobian = motionJacobian_t::Zero();

    /**
     * @brief Jacobian matrix of the process noise.
     * @details The formula used for the calculation of the process noise
     * Jacobian takes the form
     * \f[L = \begin{bmatrix}\Delta t \cos(\theta) & 0 \\ \Delta t \sin(\theta)
     * & 0 \\ 0 & \Delta t \end{bmatrix}, \f] where \f$\Delta t\f$ denotes the
     * sample period; and \f$\theta\f$ denotes the heading (orientation) of the
     * ego robot. measurement_noise. See EKF::prediction for information on the
     * motion model from which this was derived.
     */
    processJacobian_t process_jacobian = processJacobian_t::Zero();

    /**
     * @brief Jacobian of the measurement model: 2 x 5 matrix.
     * @details The formula used for the calculation of the Jacobian of the
     * measurement matrix between ego vehicle \f$i\f$ and measured vehicle
     * \f$j\f$ take the form
     * \f[ H = \begin{bmatrix} \frac{-\Delta x}{d} & \frac{-\Delta y}{d} & 0 &
     * \frac{\Delta x}{d} & \frac{\Delta y}{d} \\ \frac{\Delta y}{d^2} &
     * \frac{-\Delta x}{d^2} & -1 & \frac{-\Delta y}{d^2} & \frac{\Delta x}{d^2}
     * \end{bmatrix} \f] where \f$\Delta x = x_j - x_i\f$; \f$\Delta y = y_j
     * - y_i\f$; and \f$\Delta d = \sqrt{\Delta x^2 + \Delta y^2}\f$.
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
    covariance_t precision_matrix = error_covariance.inverse();
  };

  /**
   * @brief Houses all estimation parameters for all robots.
   */
  std::vector<EstimationParameters> robot_parameters;

  /**
   * @brief Houses all estimation parameters for all landmarks.
   */
  std::vector<EstimationParameters> landmark_parameters;

  void motionModel(const Robot::Odometry &, EstimationParameters &,
                   const double);

  void calculateMotionJacobian(const Robot::Odometry &, EstimationParameters &,
                               const double);

  void calculateProcessJacobian(EstimationParameters &, const double);

  measurement_t measurementModel(EstimationParameters &,
                                 const EstimationParameters &);

  void calculateMeasurementJacobian(EstimationParameters &,
                                    const EstimationParameters &);

  matrix3D_t marginalise(const matrix5D_t &);

  augmentedState_t createAugmentedState(const EstimationParameters &,
                                        const EstimationParameters &);

  augmentedCovariance_t createAugmentedCovariance(const EstimationParameters &,
                                                  const EstimationParameters &);

  void normaliseAngle(double &);

  huberMeasurementWeights_t
  HuberMeasurement(const measurement_t &, const huberMeasurementThresholds_t &);

  huberStateWeights_t HuberState(const augmentedState_t &,
                                 const huberStateThresholds_t &);

  measurement_t
  calculateNormalisedMeasurementResidual(const EstimationParameters &);

  augmentedState_t
  calculateNormalisedEstimationResidual(const EstimationParameters &);

public:
  explicit Filter(DataHandler &data);
  virtual ~Filter();
};

#endif // INCLUDE_SRC_FILTER_H_
