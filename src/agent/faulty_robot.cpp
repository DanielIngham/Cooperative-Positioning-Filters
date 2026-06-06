/**
 * @file faulty_robot.hpp
 * @author Daniel Ingham
 * @date 2026-06-06
 */
#include "CL/agent/faulty_robot.hpp"
#include "CL/agent/robot.hpp"
#include "CL/common/types.hpp"
#include "CL/utils/utils.hpp"

namespace CL {
FaultyRobot::FaultyRobot(const utias::mrclam::Robot &data) : Robot(data) {
  EstimationParameters &prior{estimates_.front()};

  prior.state_estimate(X) += 5;
  prior.state_estimate(Y) += 5;
  prior.state_estimate(ORIENTATION) += .5;
  CL::utils::normaliseAngle(prior.state_estimate(ORIENTATION));
}

} // namespace CL
