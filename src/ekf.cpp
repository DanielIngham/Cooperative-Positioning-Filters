#include "ekf.h"

/**
 * @brief EKF class constructor.
 * @param[in] data Class containing all robot data.
 * @details The data contained in the class includes:
 * - All robots states, odometry readings, and measurements.
 * - All landmarks positions.
 * - All robot sensors errors statistics.
 */
EKF::EKF(DataHandler &data) : data_(data) {}

/**
 * @brief Default destructor.
 */
EKF::~EKF() {}

/**
 * @brief Performs robot state inference using the EKF bayesian inference
 * framework for all robots provided.
 */
void EKF::peformInference() {}
