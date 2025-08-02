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

  /* Check for SAVE_INPUT definition that determines if the input data should
   * be plot.*/
#ifdef SAVE_INPUT
  data.saveExtractedData();
  data.plotExtractedData();
#endif // SAVE_INPUT

  /* The Inference plots are always saved.  */
  data.saveInferenceData();
  data.plotInferenceData();

  delete filter;

  return 0;
}
