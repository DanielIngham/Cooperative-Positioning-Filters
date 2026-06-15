/**
 * @file meas_data.hpp
 * @author Daniel Ingham
 * @date 2026-06-14
 */
#pragma once

#include <cstddef>
/**
 * Data structure housing the range and bearing measurements taken by mobile
 * agents in the VANET.
 */
class MeasData {
public:
  MeasData() = delete;
  MeasData(MeasData &&) = default;
  MeasData(const MeasData &) = default;
  MeasData &operator=(MeasData &&) = default;
  MeasData &operator=(const MeasData &) = default;
  ~MeasData() = default;

  /**
   * Constructor.
   * @param time Timestamp of measurement.
   * @param range Range in meters to observation.
   * @param bearing Bearing in radians of the observation from the vehicles
   * heading.
   * @param barcode Barcode of the observed agent.
   */
  MeasData(double time, double range, double bearing, size_t barcode);

  /**
   * Get the measurement timestamp.
   * @returns A timestamp.
   */
  double time() const;
  /**
   * Get the range measurement.
   * @returns A range measurement in meters.
   */
  double range() const;
  /**
   * Get the bearing measurement.
   * @returns A bearing measurement in radians.
   */
  double bearing() const;
  /**
   * Get the observed agents barcode.
   * @returns The barcode of the agent observed.
   */
  size_t barcode() const;

private:
  /** Timestamp of the measurement. */
  double time_{};
  /** Range in meters to the observed agent. */
  double range_{};
  /** Bearing in radians of the observation from the vehicles heading. */
  double bearing_{};
  /** Barcode of the observed agent. */
  size_t barcode_{};
};
