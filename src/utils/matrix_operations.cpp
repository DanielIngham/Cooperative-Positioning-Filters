/**
 * @file matrix_operations.cpp
 */
#include "CL/utils/matrix_operations.hpp"
#include "CL/common/estimation_parameters.hpp"
#include <iostream>

namespace CL {
matrix3D_t MatrixOperations::marginalise(const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight{
      matrix_6d.bottomRightCorner<total_states, total_states>()};

  matrix3D_t bottomRightMatrixInverse{computePseudoInverse(bottomRight)};

  matrix3D_t matrix_3d{
      matrix_6d.topLeftCorner<total_states, total_states>() -
      matrix_6d.topRightCorner<total_states, total_states>() *
          bottomRightMatrixInverse *
          matrix_6d.bottomLeftCorner<total_states, total_states>()};

  return matrix_3d;
}

state_t MatrixOperations::marginalise(const vector6D_t &vector_6d,
                                      const matrix6D_t &matrix_6d) {

  matrix3D_t bottomRight{
      matrix_6d.bottomRightCorner<total_states, total_states>()};

  matrix3D_t bottomRightInverse{computePseudoInverse(bottomRight)};

  state_t new_vector{vector_6d.head<total_states>() -
                     matrix_6d.topRightCorner<total_states, total_states>() *
                         bottomRightInverse * vector_6d.tail<total_states>()};

  return new_vector;
}

matrix3D_t MatrixOperations::computePseudoInverse(const matrix3D_t &matrix_3d) {

  /* Compute pseudo-inverse using SVD */
  Eigen::JacobiSVD<matrix3D_t> svd{matrix_3d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV};

  /* Set tolerance for singular values (adjust as needed) */
  static constexpr double tolerance{1e-10};

  /* Get singular values and compute pseudo-inverse */
  auto singular_values{svd.singularValues()};
  Eigen::VectorXd singular_values_inv{singular_values.size()};

  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  /* Reconstruct pseudo-inverse */
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}

matrix6D_t MatrixOperations::computePseudoInverse(const matrix6D_t &matrix_6d) {

  /* Compute pseudo-inverse using SVD */
  Eigen::JacobiSVD<matrix6D_t> svd{matrix_6d,
                                   Eigen::ComputeFullU | Eigen::ComputeFullV};

  /* Set tolerance for singular values (adjust as needed) */
  static constexpr double tolerance{1e-10};

  /* Get singular values and compute pseudo-inverse */
  auto singular_values{svd.singularValues()};
  Eigen::VectorXd singular_values_inv{singular_values.size()};

  for (int i{}; i < singular_values.size(); ++i) {
    if (singular_values(i) > tolerance) {
      singular_values_inv(i) = 1.0 / singular_values(i);
    } else {
      singular_values_inv(i) = 0.0;
    }
  }

  /* Reconstruct pseudo-inverse. */
  return svd.matrixV() * singular_values_inv.asDiagonal() *
         svd.matrixU().transpose();
}

augmentedInformation_t
MatrixOperations::createAugmentedVector(const state_t &ego_robot,
                                        const state_t &other_agent) {

  augmentedInformation_t augmented_vector{augmentedState_t::Zero()};

  augmented_vector.head<total_states>() = ego_robot;
  augmented_vector.tail<total_states>() = other_agent;

  return augmented_vector;
}

augmentedCovariance_t
MatrixOperations::createAugmentedMatrix(const covariance_t &ego_robot,
                                        const covariance_t &other_agent) {

  augmentedCovariance_t matrix{augmentedCovariance_t::Zero()};

  matrix.topLeftCorner<total_states, total_states>() = ego_robot;

  matrix.bottomRightCorner<total_states, total_states>() = other_agent;

  return matrix;
}

measurement_t MatrixOperations::normaliseInnovation(
    const measurement_t &innovation,
    const measurementCovariance_t &covariance) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<measurementCovariance_t> llt{covariance};

  if (llt.info() != Eigen::Success) {
    std::cerr << covariance << std::endl;
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  measurementCovariance_t covariance_sqrt{llt.matrixL()};

  measurement_t normalised_innovation{covariance_sqrt.inverse() * innovation};

  return normalised_innovation;
}

measurement_t MatrixOperations::unnormaliseInnovation(
    const measurement_t &normalised_innovation,
    const measurementCovariance_t &covariance) {

  /* Calculate the Cholesky of the innovation */
  Eigen::LLT<measurementCovariance_t> llt{covariance};

  if (llt.info() != Eigen::Success) {
    std::cout << covariance << std::endl;
    throw std::runtime_error(
        "An error has occurred with calculating the Cholesky decomposition of "
        "the innovation error covariance");
  }

  measurementCovariance_t covariance_sqrt{llt.matrixL()};

  measurement_t unnormalised_innovation{covariance_sqrt *
                                        normalised_innovation};

  return unnormalised_innovation;
}

augmentedState_t MatrixOperations::calculateNormalisedEstimationResidual(
    const EstimationParameters &filter) {

  /* Calculate the mean of the estimation residual (innovation). */
  static constexpr double regularisation{1e-3};

  Eigen::LLT<augmentedCovariance_t> error_covariance_cholesky(
      filter.kalman_gain * filter.innovation_covariance *
          filter.kalman_gain.transpose() +
      augmentedCovariance_t::Identity() * regularisation);

  if (error_covariance_cholesky.info() != Eigen::Success) {

    throw std::runtime_error("[4] An error has occurred with calculating the "
                             "Cholesky decomposition of "
                             "the estimation error covariance");
  }

  augmentedCovariance_t error_covariance_cholesky_matrix{
      error_covariance_cholesky.matrixL()};

  augmentedState_t normalised_error_residual{
      error_covariance_cholesky_matrix.inverse() * filter.estimation_residual};

  return normalised_error_residual;
}
} // namespace CL
