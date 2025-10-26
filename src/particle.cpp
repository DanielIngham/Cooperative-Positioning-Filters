#include "particle.h"
#include "types.h"
#include <DataHandler.h>

namespace Filters {

Particle::Particle(size_t samples, Data::Handler &data) : Filter{data} {
  std::random_device rd;
  gen_.seed(rd());

  /* TODO: Populate the samples with the prior. */
  for (const auto &[id, parameters] : robot_parameters) {
    auto result{particles_.emplace(
        id, Particles{samples, parameters.front().state_estimate})};

    if (!result.second) {
      throw std::runtime_error("Unable to add Robot with ID " + id +
                               " to robot parameters");
    }
  }
}

void Particle::prediction(const Data::Robot::Odometry &odometry,
                          EstimationParameters &ego) {

  const double sample_period{data_.getSamplePeriod()};

  Particles &particles{particles_.at(ego.id)};
  particles.propagate(odometry, ego, sample_period, gen_);
}

void Particle::correction(EstimationParameters &ego,
                          const EstimationParameters &agent) {}

/**
 * Initialises all the particles to the known prior state and sets the weights
 * equally.
 * @note In this implementation, it is assumed that the prior state is known,
 * and therefore the initial samples are all set to the known prior instead of
 * being drawn from a prior distribution.
 */
Particle::Particles::Particles(const size_t samples, const state_t &prior) {
  const double inital_weight{1. / samples};
  for (size_t i{}; i < samples; ++i) {
    samples_.emplace_back(prior, inital_weight);
  }
}

void Particle::Particles::propagate(const Data::Robot::Odometry &odometry,
                                    const EstimationParameters &ego,
                                    const double sample_period,
                                    std::mt19937 &gen) {

  for (auto &[state, weight] : samples_) {
    /* Update sample. */
    const input_t noise{
        sampleMultivariateNormal(input_t::Zero(), ego.process_noise, gen)};

    motionModel(odometry, state, noise, sample_period);
  }
}

Eigen::VectorXd
Particle::Particles::sampleMultivariateNormal(const Eigen::VectorXd &mean,
                                              const Eigen::MatrixXd &cov,
                                              std::mt19937 &gen) {

  const size_t dim{total_inputs};

  Eigen::LLT<Eigen::MatrixXd> llt{cov};
  Eigen::MatrixXd L{llt.matrixL()};

  std::normal_distribution<double> standard_normal{0.0, 1.0};
  Eigen::VectorXd z{dim};

  for (int i{}; i < dim; ++i)
    z(i) = standard_normal(gen);

  return mean + L * z;
}

void Particle::Particles::motionModel(const Data::Robot::Odometry &odometry,
                                      state_t &state, const input_t &noise,
                                      const double sample_period) {

  state << state(X) + (odometry.forward_velocity + noise(FORWARD_VELOCITY)) *
                          sample_period * std::cos(state(ORIENTATION)),
      state(Y) + (odometry.forward_velocity + noise(FORWARD_VELOCITY)) *
                     sample_period * std::sin(state(ORIENTATION)),
      state(ORIENTATION) +
          (odometry.angular_velocity + noise(ANGULAR_VELOCITY)) * sample_period;
}

} // namespace Filters
