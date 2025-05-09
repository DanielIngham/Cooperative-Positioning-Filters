#ifndef INCLUDE_INCLUDE_GBP_H_
#define INCLUDE_INCLUDE_GBP_H_

#include "filter.h"

#include <DataHandler/DataHandler.h>
/**
 * @class GBP
 * @brief Gaussian Belief Propagation
 */
class GBP : public Filter {
public:
  GBP(DataHandler &);
  ~GBP();

private:
};
#endif // INCLUDE_INCLUDE_GBP_H_
