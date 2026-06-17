/**
 * @file meas_data.cpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */
#include "CL/sensors/meas_data.hpp"

namespace CL::sensors {

MeasData::MeasData(double time, double range, double bearing, size_t barcode)
    : time_{time}, barcode_{barcode}, vec_{range, bearing} {}

double MeasData::time() const { return time_; }

measurement_t MeasData::vec() const { return vec_; }

size_t MeasData::barcode() const { return barcode_; }
} // namespace CL::sensors
