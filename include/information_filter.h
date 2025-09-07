#ifndef INCLUDE_INCLUDE_INFORMATION_H_
#define INCLUDE_INCLUDE_INFORMATION_H_

#include "filter.h"

#include <DataHandler.h>

namespace Filter {
class InformationFilter : public Filter {
public:
  InformationFilter(Data::Handler &data);
  ~InformationFilter() override;

  void prediction(const Data::Robot::Odometry &,
                  EstimationParameters &) override;

  void correction(EstimationParameters &,
                  const EstimationParameters &) override;

private:
};
} // namespace Filter
#endif // INCLUDE_INCLUDE_INFORMATION_H_
