#include <ArgumentHandler.h>
#include <DataHandler.h>
#include <Plotter.h>
#include <Robot.h>

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "filter.h"

#ifdef EKF_TARGET
#include "ekf.h"
#endif // EKF_TARGET

#ifdef IEKF_TARGET
#include "iekf.h"
#endif // IEKF_TARGET

#ifdef INFO_TARGET
#include "information_filter.h"
#endif // INFO_TARGET

#ifdef CMEKF_TARGET
#include "cmekf.h"
#endif // INFO_TARGET

#ifdef PARTICLE_TARGET
#include "particle.h"
#endif

int main(int argc, char *argv[]) {

  std::cout << "Filter performing correction in ";

#ifdef COUPLED
  std::cout << "coupled ";

#elif defined(DECOUPLED)
  std::cout << "decoupled ";

#endif // DECOUPLED

#ifdef ROBUST
  std::cout << "robust " << std::endl;
#endif // ROBUST

  std::cout << "mode. " << std::endl;

  Data::Handler data;

  ArgumentHandler::setArguments(argc, argv, data);

  std::unique_ptr<Filters::Filter> filter;

#ifdef EKF_TARGET
  filter = std::make_unique<Filters::EKF>(data);

#elif defined(IEKF_TARGET)
  filter = std::make_unique<Filters::IEKF>(data);

#elif defined(INFO_TARGET)

  filter = std::make_unique<Filters::InformationFilter>(data);

#elif defined(CMEKF_TARGET)

  filter = std::make_unique<Filters::CMEKF>(data);

#elif defined(PARTICLE_TARGET)
  filter = std::make_unique<Filters::Particle>(1000, data);
#else
  throw std::runtime_error("Filter target not selected");

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

  plotter.plotPoses(robots, {Data::Type::ABSOLUTE_ERROR});

  std::cout << data.getAverageRMSE() << std::endl;
  return EXIT_SUCCESS;
}
