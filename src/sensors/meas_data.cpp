/**
 * @file meas_data.cpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */
#include "CL/sensors/meas_data.hpp"

MeasData::MeasData(double time, double range, double bearing, size_t barcode)
    : time_{time}, range_{range}, bearing_{bearing}, barcode_{barcode} {}

double MeasData::time() const { return time_; }
double MeasData::range() const { return range_; }
double MeasData::bearing() const { return bearing_; }
size_t MeasData::barcode() const { return barcode_; }
