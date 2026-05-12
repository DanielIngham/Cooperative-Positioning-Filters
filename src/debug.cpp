#include <UtiasMrclam/Plotter.hpp>

#include <UtiasMrclam/DataHandler.hpp>

#include <cctype>
#include <cstdlib>

#include "CL/inference.hpp"

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

int main(int argc, char *argv[]) {

  Data::Handler data;

  data.setDataSet(Data::Set[0]);
  // data.setSimulation(1000);

  // Filter::EKF ekf{data};
  // ekf.performInference();
  CL::Inference<CL::filter::EKF> test{data};
  test.compute();
  CL::utils::PerformanceEvaluator::populateSyncedStates(test.getRobots(), data);
  data.calculateStateError();

  Data::Plotter plot{};
  const auto &robots{data.getRobots()};
  plot.plotPoses(robots, {Data::Type::ABSOLUTE_ERROR});

  // const Data::Robot::List &robots{data.getRobots()};
  // const Data::Landmark::List &landmarks{data.getLandmarks()};

  // plot.plotPoses(robots, {Data::Type::GROUNDTRUTH});
  // plot.plotMeasurements(robots);
  // plot.plotMeasurementsVector(robots, landmarks);
  // plot.plotMeasurements(robots, {Data::Type::ERROR});
  // plot.plotMeasurements(robots);
  // plot.plotMeasurementPDFs(robots);
  // plot.plotOdometryPDFs(robots);
  return EXIT_SUCCESS;
}
