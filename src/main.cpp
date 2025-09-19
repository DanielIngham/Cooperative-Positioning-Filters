#include "Plotter.h"
#include "Robot.h"

#include <ArgumentHandler.h>
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

  Data::Plotter plotter;

#ifdef SAVE_INPUT
  data.saveExtractedData();
#endif // SAVE_INPUT

  // plotter.plotPoses({Data::Plotter::GROUNDTRUTH, Data::Plotter::SYNCED}, 1);
  // plotter.plotPoses({Data::Plotter::ABSOLUTE_ERROR});
  const auto &robots{data.getRobots()};

  plotter.plotPoses(robots);

  std::cout << data.getAverageRMSE().x << std::endl
            << data.getAverageRMSE().y << std::endl
            << data.getAverageRMSE().orientation << std::endl;
  return EXIT_SUCCESS;
}
