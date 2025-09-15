#include "Plotter.h"

#include <DataHandler.h>

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <memory>

#ifdef EKF_TARGET
#include "ekf.h"
#endif // EKF_TARGET

#ifdef IEKF_TARGET
#include "iekf.h"
#endif // IEKF_TARGET

#ifdef INFO_TARGET
#include "information_filter.h"
#endif // INFO_TARGET

int main(int argc, char *argv[]) {

  Data::Handler data;

  data.setDataSet(Data::Set[0]);

  Filter::EKF ekf{data};
  ekf.performInference();

  Data::Plotter plot{data};
  plot.plotPoses({Data::Plotter::ABSOLUTE_ERROR});
  return EXIT_SUCCESS;
}
