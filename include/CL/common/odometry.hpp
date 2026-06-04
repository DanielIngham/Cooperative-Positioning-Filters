/**
 * @file odometry.hpp
 */
#pragma once

#include <Eigen/Dense>

class Odometry {
public:
  Odometry() = default;
  Odometry(Odometry &&) = default;
  Odometry(const Odometry &) = default;
  Odometry &operator=(Odometry &&) = default;
  Odometry &operator=(const Odometry &) = default;
  ~Odometry() = default;

  typedef Eigen::Vector2d input;

private:
};
