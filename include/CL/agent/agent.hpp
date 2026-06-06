/**
 * @file agent.hpp
 */
#pragma once
#include "CL/common/estimation_parameters.hpp"

#include <UtiasMrclam/agents/Agent.hpp>
#include <UtiasMrclam/agents/Landmark.hpp>
#include <UtiasMrclam/agents/Robot.hpp>

namespace CL {
/**
 * Autonomous agent in the VANET. Agents in the VANET communicate their
 * estimates with one another over the VANET.
 */
class Agent {
public:
  Agent() = default;
  Agent(Agent &&) = default;
  Agent(const Agent &) = default;
  Agent &operator=(Agent &&) = default;
  Agent &operator=(const Agent &) = default;
  ~Agent() = default;

  /**
   * Constructor
   * @param barcode Uniquely assigned barcode.
   */
  Agent(const utias::mrclam::Agent::Barcode &data);

  /**
   * @returns The barcode of the agent.
   */
  const utias::mrclam::Agent::Barcode &getBarcode() const;

  /**
   * Request for the agents current state estimate for the provided time index.
   * @param index Current timestamp.
   * @returns the estimate of the agent for the given timestamp.
   */
  virtual const EstimationParameters &broadcastEstimate(size_t index) = 0;

protected:
  /**
   * Current state estimate of the agent.
   */
  EstimationParameters estimation_;

private:
  /**
   * Agents unique barcode.
   */
  utias::mrclam::Agent::Barcode barcode_;
};
} // namespace CL
