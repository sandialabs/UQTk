#include <gtest/gtest.h>

#include "RosenbrockFunction.h"

#include "MUQ/Optimization/ModPieceCostFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(ModPieceCostFunctionTests, RosenbrockCost) {
  // the Rosenbrock cost function
  auto rosen = std::make_shared<RosenbrockModPiece>();
  auto cost = std::make_shared<ModPieceCostFunction>(rosen);

  // choose a random point
  const Eigen::VectorXd x = Eigen::Vector2d::Random();
  const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 100.0);

  // the true value
  const double cst = (1.0-x(0))*(1.0-x(0))+100.0*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));

  // check the cost evaluations
  EXPECT_DOUBLE_EQ(cst, cost->Evaluate(x, a).at(0) (0));
  EXPECT_DOUBLE_EQ(cst, cost->Cost(ref_vector<Eigen::VectorXd>({std::cref(x), std::cref(a)})));
  EXPECT_DOUBLE_EQ(cst, cost->Cost(x, a));

  // the true gradient
  const Eigen::Vector2d grad_true(-400.0*(x(1)-x(0)*x(0))*x(0)-2.0*(1.0-x(0)), 200.0*(x(1)-x(0)*x(0)));

  // compute the gradient
  const Eigen::VectorXd& grad_test0 = cost->Gradient(0, x, a, (Eigen::VectorXd)Eigen::VectorXd::Ones(2));
  const Eigen::VectorXd& grad_test1 = cost->Gradient(0, ref_vector<Eigen::VectorXd>({std::cref(x), std::cref(a)}), (Eigen::VectorXd)Eigen::VectorXd::Ones(2));

  EXPECT_DOUBLE_EQ((grad_true-grad_test0).norm(), 0.0);
  EXPECT_DOUBLE_EQ((grad_true-grad_test1).norm(), 0.0);
}
