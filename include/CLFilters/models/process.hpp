#pragma once

namespace Filters::Models {
class Process {
public:
  Process() = default;
  Process(Process &&) = delete;
  Process(const Process &) = delete;
  Process &operator=(Process &&) = delete;
  Process &operator=(const Process &) = delete;
  ~Process() = default;

private:
};

} // namespace Filters::Models
