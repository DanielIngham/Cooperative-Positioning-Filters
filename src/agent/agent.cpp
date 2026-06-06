/**
 * @file agent.cpp
 */
#include "CL/agent/agent.hpp"

namespace CL {

Agent::Agent(const utias::mrclam::Agent::Barcode &barcode) {
  barcode_ = barcode;
}

const utias::mrclam::Agent::Barcode &Agent::getBarcode() const {
  return barcode_;
}

} // namespace CL
