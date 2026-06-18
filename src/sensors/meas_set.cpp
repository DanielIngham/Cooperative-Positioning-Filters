#include "CL/sensors/meas_set.hpp"

namespace CL::sensors {

MeasSet::MeasSet(utias::mrclam::Robot::Measurement const &measurement,
                 measurementCovariance_t const &cov)
    : time_{measurement.time}, cov_{cov} {

  /* NOTE: This operates on the assumption that the three vectors (subjects,
   range, bearing) provided by the utias dataset handler are all the same
   length. */
  for (size_t i{}; i < measurement.subjects.size(); ++i) {
    meas_set_.emplace(time_, measurement.ranges[i], measurement.bearings[i],
                      measurement.subjects[i].val(), cov);
  }
}

double MeasSet::time() const { return time_; }

size_t MeasSet::size() const { return meas_set_.size(); }

} // namespace CL::sensors
