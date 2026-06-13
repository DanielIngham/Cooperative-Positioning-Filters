#include "CL/sensors/odom_data.hpp"
#include "CL/common/types.hpp"

namespace CL::sensors {

OdomData::OdomData(double time, double vel_f, double vel_w,
                   processCovariance_t cov)
    : time_{time}, cov_{cov} {
  input_ << vel_f, vel_w;
}

double OdomData::time() const { return time_; }

input_t const &OdomData::input() const { return input_; }

processCovariance_t &OdomData::cov() { return cov_; }

} // namespace CL::sensors
