#include <gtest/gtest.h>

#include "MUQ/Modeling/CombineVectors.h"

using namespace muq::Modeling;

TEST(CombineVectorsTests, Basic) {
  auto combine = std::make_shared<CombineVectors>(Eigen::Vector2i(2, 3));

  const Eigen::VectorXd x1 = Eigen::Vector2d::Random();
  const Eigen::VectorXd x2 = Eigen::Vector3d::Random();

  const Eigen::VectorXd x = combine->Evaluate(x1, x2) [0];

  EXPECT_EQ(x.size() , 5);
  EXPECT_NEAR((x.segment(0, 2)-x1).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((x.segment(2, 3)-x2).norm(), 0.0, 1.0e-10);

  {
    const Eigen::VectorXd v = Eigen::VectorXd::Random(5);
    const Eigen::VectorXd grad = combine->Gradient(0, 0, x1, x2, v);
    const Eigen::VectorXd gradFD = combine->GradientByFD(0, 0, ref_vector<Eigen::VectorXd>({x1, x2}), v);
    EXPECT_NEAR((grad-gradFD).norm(), 0.0, 1.0e-8);
  }

  {
    const Eigen::VectorXd v = Eigen::VectorXd::Random(5);
    const Eigen::VectorXd grad = combine->Gradient(0, 1, x1, x2, v);
    const Eigen::VectorXd gradFD = combine->GradientByFD(0, 1, ref_vector<Eigen::VectorXd>({x1, x2}), v);
    EXPECT_NEAR((grad-gradFD).norm(), 0.0, 1.0e-8);
  }

  {
    const Eigen::MatrixXd jac = combine->Jacobian(0, 0, x1, x2);
    const Eigen::MatrixXd jacFD = combine->JacobianByFD(0, 0, ref_vector<Eigen::VectorXd>({x1, x2}));
    EXPECT_NEAR((jac-jacFD).norm(), 0.0, 1.0e-8);
  }

  {
    const Eigen::MatrixXd jac = combine->Jacobian(0, 1, x1, x2);
    const Eigen::MatrixXd jacFD = combine->JacobianByFD(0, 1, ref_vector<Eigen::VectorXd>({x1, x2}));
    EXPECT_NEAR((jac-jacFD).norm(), 0.0, 1.0e-8);
  }

  {
    const Eigen::VectorXd v = Eigen::VectorXd::Random(2);
    const Eigen::VectorXd app = combine->ApplyJacobian(0, 0, x1, x2, v);
    const Eigen::VectorXd appFD = combine->ApplyJacobianByFD(0, 0, ref_vector<Eigen::VectorXd>({x1, x2}), v);
    EXPECT_NEAR((app-appFD).norm(), 0.0, 1.0e-8);
  }

  {
    const Eigen::VectorXd v = Eigen::VectorXd::Random(3);
    const Eigen::VectorXd app = combine->ApplyJacobian(0, 1, x1, x2, v);
    const Eigen::VectorXd appFD = combine->ApplyJacobianByFD(0, 1, ref_vector<Eigen::VectorXd>({x1, x2}), v);
    EXPECT_NEAR((app-appFD).norm(), 0.0, 1.0e-8);
  }
}
