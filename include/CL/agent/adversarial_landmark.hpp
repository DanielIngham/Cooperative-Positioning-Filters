/**
 * @file adversarial_landmark.hpp
 * @author Daniel Ingham
 * @date 2026-06-06
 */
#pragma once

#include "CL/agent/landmark.hpp"

namespace CL {
/**
 * Landmark agent that delibrately shares incorrect data over the VANET to
 * induce instability connected cooperative agents.
 */
class AdversarialLandmark : public Landmark {
public:
  AdversarialLandmark() = default;
  AdversarialLandmark(AdversarialLandmark &&) = default;
  AdversarialLandmark(const AdversarialLandmark &) = default;
  AdversarialLandmark &operator=(AdversarialLandmark &&) = default;
  AdversarialLandmark &operator=(const AdversarialLandmark &) = default;
  ~AdversarialLandmark() = default;

  /**
   * Augments the agents true information for malicious intent.
   * @param data landmark data set up data.
   */
  AdversarialLandmark(const utias::mrclam::Landmark &data);

  /**
   * Broadcasts falsified information regarding the agents state.
   * @param index the current time index.
   * @returns Falsified posterior state of the landmark.
   */
  const EstimationParameters &broadcastEstimate(size_t index) override;

private:
};
} // namespace CL
