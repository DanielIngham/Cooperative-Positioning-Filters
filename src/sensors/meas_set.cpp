#include "CL/sensors/meas_set.hpp"
#include "CL/sensors/meas_data.hpp"
#include <cassert>

namespace CL::sensors {

MeasSet::MeasSet(utias::mrclam::Robot::Measurement const &measurement)
    : time_{measurement.time} {

  /* NOTE: This operates on the assumption that the three vectors (subjects,
   range, bearing) provided by the utias dataset handler are all the same
   length. */
  for (size_t i{}; i < measurement.subjects.size(); ++i) {
    meas_set_.emplace(time_, measurement.ranges[i], measurement.bearings[i],
                      measurement.subjects[i].val());
  }
}

double MeasSet::time() const { return time_; }

auto MeasSet::begin() const { return meas_set_.begin(); }

auto MeasSet::end() const { return meas_set_.end(); }

} // namespace CL::sensors
