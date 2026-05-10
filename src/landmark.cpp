#include "CL/landmark.hpp"

namespace CL {

Landmark::Landmark(const Data::Landmark &data) {
  estimation_ = {.id = data.id(), .barcode = data.barcode()};

  estimation_.state_estimate << data.x(), data.y(), 0.0;

  /* The landmark only has two states: x and y coordintate.
   * NOTE: Although the landmark only has two states, the same data structure
   * is used for both the robot and landmark for compatibility with the
   * EFK::correction function, since only the x and y coordinate and thier
   * corresponding error covariances are used for the measurement update. */
  estimation_.error_covariance.diagonal().topRows(total_states - 1)
      << data.x_std_dev() * data.x_std_dev(),
      data.y_std_dev() * data.y_std_dev();

  /* NOTE: The landmarks don't have an estimate for thier orientation. So the
   * third diagonal element in the precsion matrix is never used. It is kept
   * for the generality of the implementation of the measurement update step:
   * both robot and landmarks use the same data structure, so they can be used
   * in the same measurement correction function.
   */
  estimation_.precision_matrix = estimation_.error_covariance.inverse();

  estimation_.precision_matrix(ORIENTATION, ORIENTATION) = 1;

  estimation_.information_vector =
      estimation_.precision_matrix * estimation_.state_estimate;
};

} // namespace CL
