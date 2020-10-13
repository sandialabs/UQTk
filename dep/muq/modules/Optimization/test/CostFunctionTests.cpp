#include <gtest/gtest.h>

#include "RosenbrockFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

TEST(CostFunctionTests, RosenbrockCost) {

  // the Rosenbrock cost function
  auto rosen = std::make_shared<RosenbrockFunction>();

  // choose a random point
  const Eigen::VectorXd x = Eigen::Vector2d::Random();
  const Eigen::VectorXd a = Eigen::VectorXd::Constant(1, 5.0);

  // the true value
  const double cst = (1.0-x(0))*(1.0-x(0))+a[0]*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0));

  // check the cost evaluations
  EXPECT_DOUBLE_EQ(cst, rosen->Cost(x));

  // the true gradient
  const Eigen::Vector2d grad_true(-4.0*a[0]*(x(1)-x(0)*x(0))*x(0)-2.0*(1.0-x(0)), 2.0*a[0]*(x(1)-x(0)*x(0)));

  // compute the gradient
  const Eigen::VectorXd& grad_test0 =
    rosen->Gradient(0, x, (Eigen::VectorXd)Eigen::VectorXd::Ones(1));
  
  EXPECT_DOUBLE_EQ((grad_true-grad_test0).norm(), 0.0);

  // the true hessian
  Eigen::Matrix2d hess_temp;
  hess_temp << 12.0*a[0]*x(0)*x(0)-4.0*a[0]*x(1)+2.0, -4.0*a[0]*x(0),
               -4.0*a[0]*x(0), 2.0*a[0];
  const Eigen::Matrix2d hess_true(hess_temp);

  // compute the Hessian
  std::vector<Eigen::VectorXd> input;
  input.push_back(x);

  const Eigen::MatrixXd& hess_test0 = rosen->Hessian(0, input);

  EXPECT_NEAR((hess_true-hess_test0).norm(), 0.0, 1.0e-5);

  // Test the Hessian action
  Eigen::Vector2d vec(-12.3, 34.6);


  const Eigen::VectorXd& hessAction_true = hess_true*vec;
  const Eigen::VectorXd& hessAction_test = rosen->ApplyHessian(0, input, vec);

  
  EXPECT_NEAR((hessAction_true-hessAction_test).norm(), 0.0, 3.0e-4);
  
}

