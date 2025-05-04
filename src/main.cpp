#include "ekf.h"

#include <DataHandler/DataHandler.h>
#include <iostream>

int main() {
  std::cout << "Hello World" << std::endl;
  DataHandler data(7000U, 0.02, 5U, 15U);
  data.saveExtractedData();
  data.plotExtractedData();

  EKF ekf(data);
  ekf.peformInference();
  data.saveStateError();
  data.plotInferenceError();

  return 0;
}
