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
#include <map>

#include "Agent.h"
#include "estimation_parameters.h"
#include "types.h"

namespace Filters {
class Filter {
public:
  explicit Filter(Data::Handler &data);
  Filter(Filter &&) = default;
  Filter(const Filter &) = default;
  Filter &operator=(Filter &&) = delete;
  Filter &operator=(const Filter &) = delete;
  virtual ~Filter();

  void performInference();

  virtual void prediction(const Data::Robot::Odometry &,
                          EstimationParameters &) = 0;

  virtual void correction(EstimationParameters &,
                          const EstimationParameters &) = 0;

private:
protected:
  /**
   * @brief Data class housing all the data pertaining to cooperative
   * localisation (positioning).
   */
  Data::Handler &data_;

  /**
   * @brief Houses all estimation parameters for all robots.
   */
  std::map<Data::Agent::ID, EstimationParameters> robot_parameters;

  /**
   * @brief Houses all estimation parameters for all landmarks.
   */
  std::map<Data::Agent::ID, EstimationParameters> landmark_parameters;

  /**
   * @brief The thresholds for the huber measurement cost function.
   * @details The vector is of the for: range and bearing.
   */
  const huberMeasurementThresholds_t measurement_thresholds{
      (huberMeasurementThresholds_t() << 0.2, 0.01).finished()};

  /**
   * @brief The thresholds for the huber state cost function.
   * @details The vector is of the for: ego x, ego y, ego orientation, agent x,
   * agent y.
   */
  const huberStateThresholds_t state_thresholds{
      (huberStateThresholds_t() << 0.15, 0.154, 0.255, 0.0104, 0.0104, 0.0)
          .finished()};

  static void motionModel(const Data::Robot::Odometry &, EstimationParameters &,
                          const double);

  static void calculateMotionJacobian(const Data::Robot::Odometry &,
                                      EstimationParameters &, const double);

  static void calculateProcessJacobian(EstimationParameters &, const double);

  [[nodiscard]] static measurement_t
  measurementModel(EstimationParameters &, const EstimationParameters &);

  static void calculateMeasurementJacobian(EstimationParameters &,
                                           const EstimationParameters &);

  [[nodiscard]] static measurementJacobian_t
  egoMeasurementJacobian(const EstimationParameters &,
                         const EstimationParameters &);

  [[nodiscard]] static measurementJacobian_t
  agentMeasurementJacobian(const EstimationParameters &,
                           const EstimationParameters &);

  [[nodiscard]] static matrix3D_t marginalise(const matrix6D_t &);

  [[nodiscard]] static state_t marginalise(const vector6D_t &,
                                           const matrix6D_t &);

  [[nodiscard]] static augmentedInformation_t
  createAugmentedVector(const state_t &, const state_t &);

  [[nodiscard]] static augmentedCovariance_t
  createAugmentedMatrix(const covariance_t &, const covariance_t &);

  [[nodiscard]] static measurement_t
  calculateNormalisedInnovation(const EstimationParameters &);

  [[nodiscard]] static augmentedState_t
  calculateNormalisedEstimationResidual(const EstimationParameters &);

  huberMeasurementWeights_t
  HuberMeasurement(const measurement_t &, const huberMeasurementThresholds_t &);

  huberStateWeights_t HuberState(const augmentedState_t &,
                                 const huberStateThresholds_t &);

  [[nodiscard]] static matrix3D_t computePseudoInverse(const matrix3D_t &);
  [[nodiscard]] static matrix6D_t computePseudoInverse(const matrix6D_t &);

  [[nodiscard]] EstimationParameters const *
  getEstimationParameters(const Data::Agent::Barcode &barcode) const;

  virtual void processMeasurements(Data::Robot::List &robots, size_t index);
};
} // namespace Filters

#endif // INCLUDE_SRC_FILTER_H_
