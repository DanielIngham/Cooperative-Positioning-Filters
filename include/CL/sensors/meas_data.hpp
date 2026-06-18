/**
 * @file meas_data.hpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */
#pragma once

#include "CL/common/types.hpp"
#include <Eigen/Dense>
#include <cstddef>

/**
 * Data structure housing the range and bearing measurements taken by mobile
 * agents in the VANET.
 */
namespace CL::sensors {

class MeasData {
public:
  MeasData() = delete;
  MeasData(MeasData &&) = default;
  MeasData(const MeasData &) = default;
  MeasData &operator=(MeasData &&) = delete;
  MeasData &operator=(const MeasData &) = delete;
  ~MeasData() = default;

  /**
   * Constructor.
   * @param time Timestamp of measurement.
   * @param range Range in meters to observation.
   * @param bearing Bearing in radians of the observation from the vehicles
   * heading.
   * @param barcode Barcode of the observed agent.
   */
  MeasData(double time, double range, double bearing, size_t barcode,
           measurementCovariance_t const &cov);

  /**
   * Get the measurement timestamp.
   * @returns A timestamp.
   */
  double time() const;

  /**
   * Get the range and bearing measurement data.
   * @returns A 2D Eigen measurement vector.
   */
  measurement_t vec() const;

  /**
   * Get the error covariance matrix of the range and bearing measurement.
   * @returns a constant reference to an eigen matrix containing the range and
   * bearing covariance matrix.
   */
  measurementCovariance_t const &cov() const;

  /**
   * Get the observed agents barcode.
   * @returns The barcode of the agent observed.
   */
  size_t barcode() const;

  bool operator<(MeasData const &rhs) const { return barcode_ < rhs.barcode_; }

private:
  /** Timestamp of the measurement. */
  double time_{};

  /** Barcode of the observed agent. */
  size_t barcode_{};

  /** Contains Range in meters to the observed agent and bearing in radians of
   * the observation from the vehicles heading. */
  measurement_t vec_{};

  measurementCovariance_t const &cov_;
};
} // namespace CL::sensors
