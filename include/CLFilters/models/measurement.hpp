#pragma once

namespace Filters::Models {
class Measurement {
public:
  Measurement() = default;
  Measurement(Measurement &&) = delete;
  Measurement(const Measurement &) = delete;
  Measurement &operator=(Measurement &&) = delete;
  Measurement &operator=(const Measurement &) = delete;
  ~Measurement() = default;

private:
};
} // namespace Filters::Models
