/**
 * @file Landmark.hpp
 */

#pragma once

#include "CL/common/estimation_parameters.hpp"
#include "UtiasMrclam/agents/Landmark.hpp"

namespace CL {

class Landmark {
public:
  Landmark() = default;
  Landmark(Landmark &&) = default;
  Landmark(const Landmark &) = default;
  Landmark &operator=(Landmark &&) = default;
  Landmark &operator=(const Landmark &) = default;
  ~Landmark() = default;

  Landmark(const Data::Landmark &data);

private:
  EstimationParameters estimation_;
};
} // namespace CL
