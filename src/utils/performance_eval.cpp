#include "CL/utils/performance_eval.hpp"

namespace CL::utils {
#if 0
0void PerformanceEvaluator::writeInnovation() {
  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};
    if (!success)
      throw std::runtime_error("Unable to create directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {
    const std::string filepath{directory + "/Robot_" + id + "_Innovation.dat"};

    /* Open the file in append mode */
    std::ofstream outFile(filepath, std::ios::app);

    if (!outFile)
      throw std::runtime_error("Failed to open file: " + filepath);

    const Data::Robot *robot{&data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {

      if (!utias::mrclam::utils::getMeasurement(robot, k))
        continue;

      outFile << parameter_list.at(k).innovation[RANGE] << '\t'
              << parameter_list.at(k).innovation[BEARING] << '\n';
    }

    outFile.close();
  }
}

void Filter::writeNormalisedInnovation() {

  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};

    if (!success)
      throw std::runtime_error("Unable to open directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {

    const std::string filepath{directory + "/Robot_" + id +
                               "_Normalised_Innovation.dat"};

    const Data::Robot *robot{&data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {
      if (!utias::mrclam::utils::getMeasurement(robot, k))
        continue;

      const measurement_t normalised_innovation{
          normaliseInnovation(parameter_list.at(k).innovation,
                              parameter_list.at(k).innovation_covariance)};

      /* Open the file in append mode */
      std::ofstream outFile(filepath, std::ios::app);
      if (!outFile)
        throw std::runtime_error("Failed to open file: " + filepath);

      outFile << normalised_innovation[RANGE] << '\t'
              << normalised_innovation[BEARING] << '\n';

      outFile.close();
    }
  }
}

void PerformanceEvaluator::writeNEES() {

  const std::string directory{data_.getDataInferenceDirectory()};

  if (!std::filesystem::exists(directory)) {
    bool success{std::filesystem::create_directories(directory)};
    if (!success)
      throw std::runtime_error("Unable to create directory: " + directory);
  }

  for (const auto &[id, parameter_list] : robot_parameters) {
    const std::string file_path{directory + "/Robot_" + id + "_NEES.dat"};

    const Data::Robot &robot{data_.getRobot(id)};

    for (size_t k{1U}; k < data_.getNumberOfSyncedDatapoints(); k++) {

      if (!utias::mrclam::utils::getMeasurement(&robot, k))
        continue;

      const EstimationParameters &parameters{parameter_list.at(k)};

      state_t groundtruth;
      groundtruth << robot.groundtruth.states.at(k).x,
          robot.groundtruth.states.at(k).y,
          robot.groundtruth.states.at(k).orientation;

      state_t error{groundtruth - parameters.state_estimate};
      utias::mrclam::utils::normaliseAngle(error(ORIENTATION));

      const double nees{error.transpose() *
                        parameters.error_covariance.inverse() * error};

      /* Open the file in append mode */
      std::ofstream outFile(file_path, std::ios::app);

      if (!outFile)
        throw std::runtime_error("Failed to open file " + file_path);

      outFile << nees << '\n';

      outFile.close();
    }
  }
}
#endif // 0
} // namespace CL::utils
