#pragma once

#include "CL/agent/landmark.hpp"

namespace CL {
class AdversarialLandmark : public Landmark {
public:
  AdversarialLandmark() = default;
  AdversarialLandmark(AdversarialLandmark &&) = default;
  AdversarialLandmark(const AdversarialLandmark &) = default;
  AdversarialLandmark &operator=(AdversarialLandmark &&) = default;
  AdversarialLandmark &operator=(const AdversarialLandmark &) = default;
  ~AdversarialLandmark() = default;

  AdversarialLandmark(const Data::Landmark &data);

  const EstimationParameters &broadcastEstimate(size_t index) override;

private:
};
} // namespace CL
