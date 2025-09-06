#include "Plotter.h"
#include <memory>
#ifdef EKF_TARGET
#include "ekf.h"
#endif // EKF

#ifdef IEKF_TARGET
#include "iekf.h"
#endif // IEKF

#ifdef INFO_TARGET
#include "information_filter.h"
#endif // INFO

#include <ArgumentHandler.h>
#include <DataHandler.h>
#include <cctype>
#include <iostream>

using std::make_unique;

int main(int argc, char *argv[]) {

  std::cout << "Filter performing correction in ";

#ifdef COUPLED
  std::cout << "coupled ";

#elif defined(DECOUPLED)
  std::cout << "decoupled ";

#endif // DECOUPLED

#ifdef ROBUST
  std::cout << " robust " << std::endl;
#endif // ROBUST

  std::cout << " mode. " << std::endl;

  Data::Handler data;

  ArgumentHandler::setArguments(argc, argv, data);

  std::unique_ptr<Filter::Filter> filter;

#ifdef EKF_TARGET
  filter = make_unique<Filter::EKF>(data);

#elif defined(IEKF_TARGET)
  filter = make_unique<Filter::IEKF>(data);

#else
  filter = make_unique<Filter::InformationFilter>(data);

#endif // EKF_TARGET

  filter->performInference();

  /* Check for SAVE_INPUT definition that determines if the input data should
   * be plot.*/

  Data::Plotter plotter(data);

#ifdef SAVE_INPUT
  data.saveExtractedData();
#endif // SAVE_INPUT

  // plotter.plotPoses({Data::Plotter::GROUNDTRUTH, Data::Plotter::SYNCED}, 1);
  plotter.plotPoses({Data::Plotter::ABSOLUTE_ERROR}, 1);

  return 0;
}
