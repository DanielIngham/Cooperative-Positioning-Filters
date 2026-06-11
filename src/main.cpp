#include <UtiasMrclam/DataHandler.hpp>
#include <UtiasMrclam/Plotter.hpp>
#include <UtiasMrclam/agents/Robot.hpp>
#include <UtiasMrclam/utils/ArgumentHandler.hpp>

#include "CL/common/config.hpp"
#include "CL/inference.hpp"

#include <cctype>
#include <cstdlib>
#include <iostream>

#include "CL/utils/performance_eval.hpp"

#ifdef EKF_TARGET
#include "CL/filters/ekf.hpp"
#endif // EKF_TARGET

#ifdef IEKF_TARGET
#include "CL/filters/iekf.hpp"
#endif // IEKF_TARGET

#ifdef INFO_TARGET
#include "CL/filters/information_filter.hpp"
#endif // INFO_TARGET

#ifdef CMEKF_TARGET
#include "CL/filters/cmekf.hpp"
#endif // INFO_TARGET

#ifdef PARTICLE_TARGET
#include "CL/filters/particle.hpp"
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

  utias::mrclam::Handler data;

  utias::mrclam::ArgumentHandler::setArguments(argc, argv, data);

  CL::Config config{"../config/example.yaml"};

#ifdef EKF_TARGET
  CL::Inference<CL::filter::EKF> inference{data, config};

#elif defined(IEKF_TARGET)
  CL::Inference<CL::filter::IEKF> inference{data, config};

#elif defined(INFO_TARGET)
  std::cerr << "Using Information filter." << std::endl;
  CL::Inference<CL::filter::InformationFilter> inference{data, config};

#elif defined(CMEKF_TARGET)
  CL::Inference<CL::filter::CMEKF> inference{data, config};

#elif defined(PARTICLE_TARGET)
  CL::Inference<CL::filter::Particle> inference{data, config};

#else
  throw std::runtime_error("Filter target not selected");

#endif // EKF_TARGET

  inference.compute();
  CL::utils::PerformanceEvaluator::populateSyncedStates(inference.getRobots(),
                                                        data);
  data.calculateStateError();

  /* Check for SAVE_INPUT definition that determines if the input data
   * should be plot.*/

  utias::mrclam::Plotter plotter;

#ifdef SAVE_INPUT
  data.saveExtractedData();
#endif // SAVE_INPUT

  const auto &robots{data.getRobots()};
  const auto &landmarks{data.getLandmarks()};

  plotter.plotPoses(robots, {utias::mrclam::Type::ABSOLUTE_ERROR});

  std::cout << data.getAverageRMSE() << std::endl;

  return EXIT_SUCCESS;
}
