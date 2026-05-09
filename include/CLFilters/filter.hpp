/**
 * @file filter.h
 * @brief Header file of the parent class containing shared functionality amoung
 * cooperative localsiation filters.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#pragma once

#include "CLFilters/common/estimation_parameters.hpp"
#include "CLFilters/common/types.hpp"
#include "CLFilters/models/measurement.hpp"
#include "CLFilters/models/process.hpp"

#include <Eigen/Dense>
#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/agents/Agent.hpp>
#include <map>

namespace Filters {
class Filter {
public:
  Filter() = delete;
  explicit Filter(Data::Handler &data);
  Filter(Filter &&) = default;
  Filter(const Filter &) = default;
  Filter &operator=(Filter &&) = delete;
  Filter &operator=(const Filter &) = delete;
  virtual ~Filter();

  using ParameterList = std::vector<EstimationParameters>;
  using ParameterMap = std::map<Data::Agent::ID, ParameterList>;

  /**
   * @brief Performs robot state inference using the EKF bayesian inference
   * framework for all robots provided.
   */
  void performInference();

  virtual void prediction(const Data::Robot::Odometry &,
                          EstimationParameters &) = 0;

  virtual void correction(EstimationParameters &,
                          const EstimationParameters &) = 0;

  [[nodiscard]] static measurement_t
  normaliseInnovation(const measurement_t &, const measurementCovariance_t &);

  [[nodiscard]] static measurement_t
  unnormaliseInnovation(const measurement_t &, const measurementCovariance_t &);

  [[nodiscard]] static augmentedState_t
  calculateNormalisedEstimationResidual(const EstimationParameters &);

  [[nodiscard]] EstimationParameters const *
  getEstimationParameters(const Data::Agent::Barcode &barcode) const;

  void writeInnovation();
  void writeNormalisedInnovation();
  void writeNEES();

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
  std::map<Data::Agent::ID, ParameterList> robot_parameters;

  /**
   * @brief Houses all estimation parameters for all landmarks.
   */
  std::map<Data::Agent::ID, EstimationParameters> landmark_parameters;

  virtual void processMeasurements(Data::Robot::List &robots, size_t index);
};
} // namespace Filters
