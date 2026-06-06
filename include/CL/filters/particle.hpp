#pragma once
#include "CL/common/estimation_parameters.hpp"
#include "CL/common/types.hpp"
#include "CL/filters/filter.hpp"

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

  virtual EstimationParameters
  prediction(const utias::mrclam::Robot::Odometry &odometry,
             const EstimationParameters &ego, double sample_period) override;

  virtual void correction(EstimationParameters &ego,
                          const EstimationParameters &agent) override;

protected:
  struct Particles {
    /**
     * Initialises all the particles to the known prior state and sets the
     * weights equally.
     * @note In this implementation, it is assumed that the prior state is
     * known, and therefore the initial samples are all set to the known prior
     * instead of being drawn from a prior distribution.
     */
    Particles(const size_t, const state_t &);

    void propagate(const utias::mrclam::Robot::Odometry &,
                   const EstimationParameters &, const double, std::mt19937 &);

    /**
     * Calculates the new weights of the particles based on the likelihood of
     * the state given the measurement.
     * @returns A flag indicating whether the particles should be resampled.
     */
    bool reweight(const EstimationParameters &, const EstimationParameters &);

    state_t mmse();

    void resample(std::mt19937 &gen);

  private:
    /** sample and weight vector pair. */
    std::vector<std::pair<state_t, double>> samples_;

    Eigen::VectorXd sampleMultivariateNormal(const Eigen::VectorXd &mean,
                                             const Eigen::MatrixXd &cov,
                                             std::mt19937 &gen);

    void motionModel(const utias::mrclam::Robot::Odometry &odometry,
                     state_t &state, const input_t &noise,
                     const double sample_period);

    double Gaussian(const Eigen::VectorXd &mean,
                    const Eigen::MatrixXd &covariance);

    bool checkWeights();
  };

private:
  std::mt19937 gen_;

  Particles particles_;
};

} // namespace CL::filter
