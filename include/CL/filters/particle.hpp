#pragma once
#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include "CL/filters/filter.hpp"
#include "CL/sensors/meas_data.hpp"
#include "CL/sensors/odom_data.hpp"

#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <random>

namespace CL::filter {

class Particle : public Filter {
public:
  Particle() = delete;
  Particle(const EstimationParameters &prior, size_t samples = 1000);
  Particle(Particle &&) = default;
  Particle(const Particle &) = default;
  Particle &operator=(Particle &&) = delete;
  Particle &operator=(const Particle &) = delete;
  ~Particle() = default;

  virtual EstimationParameters prediction(sensors::OdomData const &odometry,
                                          EstimationParameters const &ego,
                                          double sample_period) override;

  virtual void correction(EstimationParameters &ego,
                          EstimationParameters const &agent,
                          sensors::MeasData const &meas) override;

protected:
  /**
   * Set of particle used in Particle filtering.
   */
  struct Particles {
    /**
     * Initialises all the particles to the known prior state and sets the
     * weights equally.
     * @note In this implementation, it is assumed that the prior state is
     * known, and therefore the initial samples are all set to the known prior
     * instead of being drawn from a prior distribution.
     */
    Particles(const size_t, const state_t &);

    /**
     * Propagtes the particles using the motion model.
     */
    void propagate(sensors::OdomData const &odometry,
                   EstimationParameters const &ego, double sample_period,
                   std::mt19937 &gen);

    /**
     * Calculates the new weights of the particles based on the likelihood of
     * the state given the measurement.
     * @param ego Estimate parameters of the ego vehicle.
     * @param agent Estimate parameters of the observed agent.
     * @param meas Measurement of the observed agent produced by the ego agent.
     * @returns A flag indicating whether the particles should be resampled.
     */
    bool reweight(const EstimationParameters &ego,
                  const EstimationParameters &agent, sensors::MeasData meas);

    state_t mmse();

    void resample(std::mt19937 &gen);

  private:
    /** sample and weight vector pair. */
    std::vector<std::pair<state_t, double>> samples_;

    Eigen::VectorXd sampleMultivariateNormal(const Eigen::VectorXd &mean,
                                             const Eigen::MatrixXd &cov,
                                             std::mt19937 &gen);

    void motionModel(sensors::OdomData const &odometry, state_t &state,
                     const input_t &noise, const double sample_period);

    double Gaussian(const Eigen::VectorXd &mean,
                    const Eigen::MatrixXd &covariance);

    bool checkWeights();
  };

private:
  std::mt19937 gen_;

  Particles particles_;
};

} // namespace CL::filter
