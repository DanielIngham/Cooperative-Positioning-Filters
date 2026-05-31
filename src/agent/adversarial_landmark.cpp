#include "CL/agent/adversarial_landmark.hpp"
#include "CL/common/types.hpp"

namespace CL {
AdversarialLandmark::AdversarialLandmark(const Data::Landmark &data)
    : Landmark{data} {
  state_t &falsified_state{estimation_.state_estimate};
  falsified_state(X) += 10;
  falsified_state(Y) += 10;
}

const EstimationParameters &
AdversarialLandmark::broadcastEstimate(size_t index) {
  return estimation_;
}
} // namespace CL
