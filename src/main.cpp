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

int main(int argc, char *argv[]) {

#ifdef COUPLED
  std::cout << "coupled ";

#elif defined(DECOUPLED)
  std::cout << "decoupled ";

#endif // DECOUPLED

#ifdef ROBUST
  std::cout << ", robust " << std::endl;
#endif // ROBUST

  std::cout << std::endl;

  /* Creating an instance of the DataHandler class whose output directory will
   * will be set by the respective filter.*/
  DataHandler data;

  ArgumentHandler::setArguments(argc, argv, data);

  Filter *filter;

#ifdef EKF_TARGET
  filter = new EKF(data);

#elif defined(IEKF_TARGET)
  filter = new IEKF(data);

#else
  filter = new InformationFilter(data);

#endif // EKF_TARGET

  filter->performInference();

  /* Check for environment variable IN that determines if the input data should
   * be plot.*/
  const char *plot_input = std::getenv("IN");
  if (plot_input != nullptr) {
    if (std::strcmp(plot_input, "1")) {
      data.saveExtractedData();
      data.plotExtractedData();
    }
  }

  /* The Inference plots are always saved.  */
  data.saveInferenceData();
  data.plotInferenceData();

  delete filter;

  return 0;
}
