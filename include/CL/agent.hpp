#pragma once
#include "CL/common/estimation_parameters.hpp"

#include <UtiasMrclam/agents/Agent.hpp>
#include <UtiasMrclam/agents/Landmark.hpp>
#include <UtiasMrclam/agents/Robot.hpp>

namespace CL {
class Agent {
public:
  Agent() = default;
  Agent(Agent &&) = default;
  Agent(const Agent &) = default;
  Agent &operator=(Agent &&) = default;
  Agent &operator=(const Agent &) = default;
  ~Agent() = default;

  Agent(const Data::Agent::Barcode &data);

  const Data::Agent::Barcode &getBarcode() const;

  virtual const EstimationParameters &broadcastEstimate(size_t index) = 0;

protected:
  EstimationParameters estimation_;

private:
  Data::Agent::Barcode barcode_;
};
} // namespace CL
