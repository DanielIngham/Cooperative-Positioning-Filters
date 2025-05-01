/**
 * @file ekf.cpp
 * @brief Implementation of the Extended Kalman Fitler implementation for
 * multirobot cooperative positioning.
 * @author Daniel Ingham
 * @date 2025-05-01
 */
#include "ekf.h"

/**
 * @brief EKF class constructor.
 * @param[in] data Class containing all robot data.
 * @details The data contained in the class includes:
 * - All robots states, odometry readings, and measurements.
 * - All landmarks positions.
 * - All robot sensors errors statistics.
 */
EKF::EKF(DataHandler &data) : data_(data) {
  // for (unsigned short id = 0; id < data_.getRobots().size(); id++) {
  //   this->robots.push_back(EstimationParameters());
  //   this->robots.back().error_covarince[];
  // }
}

/**
 * @brief Default destructor.
 */
EKF::~EKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void EKF::peformInference() {}
