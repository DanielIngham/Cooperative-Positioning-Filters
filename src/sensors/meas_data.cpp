/**
 * @file meas_data.cpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */
#include "CL/sensors/meas_data.hpp"

namespace CL::sensors {

MeasData::MeasData(double time, double range, double bearing, size_t barcode,
                   measurementCovariance_t const &cov)
    : time_{time}, barcode_{barcode}, vec_{range, bearing}, cov_{cov} {}

double MeasData::time() const { return time_; }

measurement_t MeasData::vec() const { return vec_; }

measurementCovariance_t const &MeasData::cov() const { return cov_; }

size_t MeasData::barcode() const { return barcode_; }
} // namespace CL::sensors
