#include "Plotter.h"

#include <DataHandler.h>

#include <cctype>
#include <cstdlib>

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
  // data.setSimulation(1000);

  Filter::EKF ekf{data};
  ekf.performInference();

  Data::Plotter plot{};
  const Data::Robot::List &robots{data.getRobots()};
  const Data::Landmark::List &landmarks{data.getLandmarks()};

  // plot.plotPoses(robots, {Data::Type::GROUNDTRUTH});
  // plot.plotMeasurements(robots);
  // plot.plotMeasurementsVector(robots, landmarks);
  // plot.plotMeasurements(robots, {Data::Type::ERROR});
  // plot.plotMeasurements(robots);
  plot.plotMeasurementPDFs(robots);
  // plot.plotOdometryPDFs(robots);
  return EXIT_SUCCESS;
}
