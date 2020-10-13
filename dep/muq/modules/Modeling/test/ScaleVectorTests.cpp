#include <gtest/gtest.h>

#include "MUQ/Modeling/ScaleVector.h"

using namespace muq::Modeling;

TEST(ScaleVectorTests, Basic) {
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(10);

  double const s = 3.0;

  auto scale = std::make_shared<ScaleVector>(s, vec.size());

  const Eigen::VectorXd& result = scale->Evaluate(vec) [0];
  EXPECT_NEAR((result-s*vec).norm(), 0.0, 1.0e-10);

  const Eigen::MatrixXd& jac = scale->Jacobian(0, 0, vec);
  EXPECT_NEAR((jac-s*Eigen::MatrixXd::Identity(vec.size(), vec.size())).norm(), 0.0, 1.0e-10);

  const Eigen::VectorXd sens = Eigen::VectorXd::Random(vec.size());
  const Eigen::VectorXd& grad = scale->Gradient(0, 0, vec, sens);
  EXPECT_NEAR((grad-jac.transpose()*sens).norm(), 0.0, 1.0e-10);

  const Eigen::VectorXd targ = Eigen::VectorXd::Random(vec.size());
  const Eigen::VectorXd& jacapp = scale->Gradient(0, 0, vec, targ);
  EXPECT_NEAR((grad-jac*targ).norm(), 0.0, 1.0e-10);
}
