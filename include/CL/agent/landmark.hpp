/**
 * @file Landmark.hpp
 */

#pragma once

#include "CL/agent/agent.hpp"
#include "CL/common/estimation_parameters.hpp"
#include "UtiasMrclam/agents/Landmark.hpp"

namespace CL {

class Landmark : public Agent {
public:
  Landmark() = default;
  Landmark(Landmark &&) = default;
  Landmark(const Landmark &) = default;
  Landmark &operator=(Landmark &&) = default;
  Landmark &operator=(const Landmark &) = default;
  ~Landmark() = default;

  Landmark(const Data::Landmark &data);

  const EstimationParameters &broadcastEstimate(size_t index) override;

private:
};
} // namespace CL
