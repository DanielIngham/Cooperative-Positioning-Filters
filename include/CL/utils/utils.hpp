/**
 * @file utils.hpp
 */
#include <cmath>
#include <numbers>
#pragma once

namespace CL::utils {

/**
 * Normalises a given angle to fall in the range [-pi, pi).
 * @param angle_rad unnormalised angle in radians
 * @returns an equivalent angle normalised to fall within the range [-pi,
 * pi).
 */
inline void normaliseAngle(double &angle_rad) {
  angle_rad -= 2.0 * std::numbers::pi *
               floor((angle_rad + std::numbers::pi) / (2.0 * std::numbers::pi));
}
} // namespace CL::utils
