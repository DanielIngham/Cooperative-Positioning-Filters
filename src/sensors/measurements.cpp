#include "CL/sensors/measurements.hpp"
#include "CL/sensors/meas_set.hpp"
#include <algorithm>

namespace apn {

Measurements::Measurements(
    std::vector<utias::mrclam::Robot::Measurement> const &measurements) {
  for (const auto &meas : measurements) {
    measurements_.emplace_back(meas);
  }
}

MeasSet *const Measurements::find(double time) {

  auto it{std::lower_bound(
      measurements_.begin(), measurements_.end(), time,
      [](const MeasSet &meas, double t) { return meas.time() < t; })};

  if (it == measurements_.end())
    return nullptr;

  /* Threshold for double floating point precision. */
  static constexpr double decimal_threshold{1e-5};
  for (size_t i{}; i < 2; ++i) {

    if (std::abs(it->time() - time) < decimal_threshold)
      return &(*it);

    if (it == measurements_.begin())
      break;

    --it;
  }

  return nullptr;
}

} // namespace apn
