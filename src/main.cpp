#include "ekf.h"
#include "iekf.h"
#include "information_filter.h"

#include <DataHandler/DataHandler.h>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>

int main() {
  /* Extract the filter environment variable */
  const char *filter_ptr = std::getenv("Filter");

  if (filter_ptr == nullptr) {
    std::cout << std::endl;
    std::cerr << "\033[1;31mERROR\033[0m: Filter not set.\n"
              << "Please set the filter using: "
              << "\033[2mmake run Filter=EKF.\033[0m" << std::endl;
    std::cout << std::endl;
    return 1;
  }

  /* Convert the string to upper case, to allow for variations input filter
   * name. */
  std::string filter(filter_ptr);
  std::transform(filter.begin(), filter.end(), filter.begin(), ::toupper);

  std::cout << "\033[1m" << filter << ": \033[0m";

#ifdef COUPLED
  std::cout << "coupled ";
#endif
#ifdef DECOUPLED
  std::cout << "decoupled ";
#endif // DECOUPLED

#ifdef ROBUST
  std::cout << ", robust " << std::endl;
#endif

  std::cout << std::endl;

  const char *dataset_number = std::getenv("Dataset");
  std::string dataset;

  /* Set the Default data date to dataset 1. */
  if (dataset_number == nullptr) {
    dataset = "MRCLAM_Dataset1";

  } else {
    dataset = "MRCLAM_Dataset" + std::string(dataset_number);
  }

  /* Creating an instance of the DataHandler class whose output directory will
   * will be set by the respective filter.*/
  DataHandler data;

  /* Check which type of filter was requested by the user and run that specific
   * one */
  if ("EKF" == filter) {
    data.setDataSet(dataset, "EKF");

    EKF ekf(data);

    ekf.performInference();

  } else if ("IEKF" == filter) {

    data.setDataSet(dataset, "IEKF");

    IEKF iekf(data);

    iekf.performInference();

  } else if ("INFO" == filter) {

    data.setDataSet(dataset, "Information_Filter");

    InformationFilter info(data);

    info.performInference();

  } else {
    std::cerr
        << "Filter specified does not exist. Here is a list of all filters:\n"
           "\tEKF, IEKF, INFO \n"
           "Type: make run Filter=EKF"
        << std::endl;
    return 1;
  }

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
  data.saveStateError();
  data.plotInferenceError();

  return 0;
}
