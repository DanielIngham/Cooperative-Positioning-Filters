/**
 * @file matrix_operations.hpp
 */
#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"

namespace CL {
class MatrixOperations {
public:
  MatrixOperations() = default;
  MatrixOperations(MatrixOperations &&) = delete;
  MatrixOperations(const MatrixOperations &) = delete;
  MatrixOperations &operator=(MatrixOperations &&) = delete;
  MatrixOperations &operator=(const MatrixOperations &) = delete;
  ~MatrixOperations() = default;

  /**
   * @brief Schur complement-based marginalisation that marginalises a 6x6
   * matrix into a 3x3 matrix.
   * @param[in] matrix_6d A 6x6 matrix.
   * @returns A 3x3 matrix.
   */
  [[nodiscard]] static matrix3D_t marginalise(const matrix6D_t &);

  /**
   * @brief Schur complement-based marginalisation that marginalises a 6x1
   * vector and 6x6 matrix into a 3x1 vector.
   * @param[in] vector_6d A 6x1 vector.
   * @param[in] matrix_6d A 6x6 matrix.
   * @returns A 3x1 vector.
   */
  [[nodiscard]] static state_t marginalise(const vector6D_t &,
                                           const matrix6D_t &);

  [[nodiscard]] static matrix3D_t computePseudoInverse(const matrix3D_t &);

  [[nodiscard]] static matrix6D_t computePseudoInverse(const matrix6D_t &);

  /**
   * @brief Combines the 3 states of the ego vehicle (x,y,orientation), with the
   * state of the agent measured to create a augmented state vector.
   * @param[in] ego_robot the structure containing the estimation parameters of
   * the ego vehicle.
   * @param[in] other_agent the structure containing the estimation parameters
   * of the agent measured by the ego vehicle.
   * @returns A 5x1 state vector.
   */
  [[nodiscard]] static augmentedInformation_t
  createAugmentedVector(const state_t &, const state_t &);

  /**
   * @brief Combines the covariance/precision matrix of the 3 states of the ego
   * vehicle (x,y,orientation), with the covariance/precision of the of the
   * agent measured.
   *
   * @param[in] ego_robot the structure containing the estimation parameters of
   * the ego vehicle.
   * @param[in] other_agent the structure containing the estimation parameters
   * of the agent measured by the ego vehicle.
   *
   * @details Cooperative Localisation (Positioning) involves robots that share
   * thier state and estimation error covariances / precision when one robot
   * measures the other. As a result, the estimation error covariance/precision
   * needs to be augmented from a 3x3 to a 6x6 matrix to house the error
   * covariance of both the ego vehicle (\f$i\f$) and the measured agent
   * (\f$j\f$):
   * \f[\mathbf{P} = \begin{bmatrix} \mathbf{P}_i & \mathbf{0} \\ \mathbf{0} &
   * \mathbf{P}_j \end{bmatrix}, \f] where \f$\mathbf{P}_i\f$ and
   * \f$\mathbf{P}_j\f$ are the estimation error covariance of the ego robot
   * \f$i\f$ and the observed agent \f$j\f$ respectively.
   */
  [[nodiscard]] static augmentedCovariance_t
  createAugmentedMatrix(const covariance_t &, const covariance_t &);

  /**
   * @brief Calculates the normalised residual of the value produced by the
   * sensor measurement and the prior estimate passed through the non-linear
   * measurement model.
   * @param[in] filter The estimation parameters of the filter.
   * @returns The normalised innovation.
   */
  [[nodiscard]] static measurement_t
  normaliseInnovation(const measurement_t &, const measurementCovariance_t &);

  [[nodiscard]] static measurement_t
  unnormaliseInnovation(const measurement_t &, const measurementCovariance_t &);

  /**
   * @brief Calculates the normalised residual of the initial estimate and
   * updated estimate.
   * @param[in] filter The estimation parameters of the filter.
   * @returns The normalised estimation residual.
   */
  [[nodiscard]] static augmentedState_t
  calculateNormalisedEstimationResidual(const EstimationParameters &);

private:
};
} // namespace CL
