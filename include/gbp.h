#ifndef INCLUDE_INCLUDE_GBP_H_
#define INCLUDE_INCLUDE_GBP_H_

#include "filter.h"

#include <DataHandler/DataHandler.h>
#include <vector>
/**
 * @class GBP
 * @brief Gaussian Belief Propagation
 */
class GBP : public Filter {
public:
  GBP(DataHandler &, const unsigned short int window_size = 2U);
  ~GBP();

  void performInference();

private:
  unsigned short int window_size_;

  std::vector<Filter::EstimationParameters> graph;

  void prediction(const Robot::Odometry &,
                  std::vector<Filter::EstimationParameters> &);
  void correction(std::vector<EstimationParameters> &,
                  const EstimationParameters &);
};
#endif // INCLUDE_INCLUDE_GBP_H_
