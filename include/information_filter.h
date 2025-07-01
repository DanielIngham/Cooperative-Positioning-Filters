#ifndef INCLUDE_INCLUDE_INFORMATION_H_
#define INCLUDE_INCLUDE_INFORMATION_H_

#include "filter.h"

#include <DataHandler/DataHandler.h>

class InformationFilter : public Filter {
public:
  InformationFilter(DataHandler &data);
  ~InformationFilter() override;

private:
  void prediction(const Robot::Odometry &, EstimationParameters &) override;
  void correction(EstimationParameters &,
                  const EstimationParameters &) override;
};

#endif // INCLUDE_INCLUDE_INFORMATION_H_
