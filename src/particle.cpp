#include "particle.h"
#include "types.h"
#include <DataHandler.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>

namespace Filters {

Particle::Particle(size_t samples, Data::Handler &data) : Filter{data} {
  std::random_device rd;
  gen_.seed(rd());

  /* Populate the samples with the prior. */
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

  ego.state_estimate = particles.mmse();
}

void Particle::correction(EstimationParameters &ego,
                          const EstimationParameters &agent) {

  Particles &particles{particles_.at(ego.id)};
  bool resample{particles.reweight(ego, agent)};

  if (resample)
    particles.resample(gen_);

  ego.state_estimate = particles.mmse();
}

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
    const input_t noise{
        sampleMultivariateNormal(input_t::Zero(), ego.process_noise, gen)};

    motionModel(odometry, state, noise, sample_period);
  }
}

/**
 * Calculates the new weights of the particles based on the likelihood of the
 * state given the measurement.
 * @returns A flag indicating whether the particles should be resampled.
 */
bool Particle::Particles::reweight(const EstimationParameters &ego,
                                   const EstimationParameters &agent) {
  for (auto &[state, weight] : samples_) {
    const measurement_t difference{
        ego.measurement - measurementModel(state, agent.state_estimate)};
    weight *= Gaussian(difference, ego.measurement_noise);
  }

  double total_weight{
      std::accumulate(samples_.begin(), samples_.end(), 0.0,
                      [](double acc, const std::pair<state_t, double> &sample) {
                        return acc + sample.second;
                      })};

  for (auto &[state, weight] : samples_) {
    weight /= total_weight;
  }

  double effective_samples{
      std::accumulate(samples_.begin(), samples_.end(), 0.0,
                      [](double acc, const std::pair<state_t, double> &sample) {
                        return acc + std::pow(sample.second, 2);
                      })};

  effective_samples = 1 / effective_samples;

  if (effective_samples < samples_.size() / 4.0)
    return true;

  return false;
}

void Particle::Particles::resample(std::mt19937 &gen) {
  std::vector<std::pair<state_t, double>> weight_map(samples_.size());

  std::inclusive_scan(samples_.begin(), samples_.end(), weight_map.begin(),
                      [](const std::pair<state_t, double> &acc,
                         const std::pair<state_t, double> &curr) {
                        return std::make_pair(curr.first,
                                              acc.second + curr.second);
                      });

  std::uniform_real_distribution<double> uniform{0., 1.};
  for (auto &[sample, weight] : samples_) {

    auto it{std::lower_bound(weight_map.begin(), weight_map.end(), uniform(gen),
                             [](const std::pair<state_t, double> &a, double b) {
                               return a.second < b;
                             })};

    sample = it->first;
    weight = 1.0 / samples_.size();
  }
}

state_t Particle::Particles::mmse() {
  state_t init{state_t::Zero()};

  return std::accumulate(
      samples_.begin(), samples_.end(), init,
      [](state_t acc, const std::pair<state_t, double> &pair) {
        return acc + pair.second * pair.first;
      });
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

double Particle::Particles::Gaussian(const Eigen::VectorXd &difference,
                                     const Eigen::MatrixXd &covariance) {
  /* Covariance matrix must be square. */
  assert(covariance.rows() == covariance.cols());
  /* Covariance matrix and mean vector must be correctly sized. */
  assert(covariance.cols() == difference.size());

  return std::exp(-0.5 * difference.transpose() * covariance.inverse() *
                  difference);
}

} // namespace Filters
