#include "CL/agent.hpp"

namespace CL {

Agent::Agent(const Data::Agent::Barcode &barcode) { barcode_ = barcode; }
const Data::Agent::Barcode &Agent::getBarcode() const { return barcode_; }

} // namespace CL
