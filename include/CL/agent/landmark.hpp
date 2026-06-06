/**
 * @file Landmark.hpp
 * @author Daniel Ingham
 * @date 2026-06-06
 */

#pragma once

#include "CL/agent/agent.hpp"
#include "CL/common/estimation_parameters.hpp"

#include <UtiasMrclam/agents/Landmark.hpp>

namespace CL {

/**
 * Landmark agent. Stationary agent, such infrustructure which should service as
 * an accurate localisation source within the VANET.
 */
class Landmark : public Agent {
public:
  Landmark() = default;
  Landmark(Landmark &&) = default;
  Landmark(const Landmark &) = default;
  Landmark &operator=(Landmark &&) = default;
  Landmark &operator=(const Landmark &) = default;
  ~Landmark() = default;

  /**
   * Constructor copies the landmark state information from the
   * dataset/simulation.
   */
  Landmark(const utias::mrclam::Landmark &data);

  /**
   * Broadcasts the parameters of Gaussian distribution of the vehicles
   * landmarks posterior state.
   * @param index the current time index.
   * @returns Posterior state of the landmark.
   */
  const EstimationParameters &broadcastEstimate(size_t index) override;

private:
};
} // namespace CL
