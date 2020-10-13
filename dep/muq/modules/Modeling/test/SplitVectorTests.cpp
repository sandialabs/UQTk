#include <gtest/gtest.h>

#include "MUQ/Modeling/SplitVector.h"

using namespace muq::Modeling;

TEST(SplitVectorTests, Basic) {
  // create a vector we want to split
  Eigen::VectorXd vec = Eigen::VectorXd::Random(10);

  // the starting index
  Eigen::VectorXi ind = Eigen::Vector3i(0, 4, 8);

  // the size of each section
  Eigen::VectorXi size = Eigen::Vector3i(4, 4, 2);

  // create the splitter
  auto split = std::make_shared<SplitVector>(ind, size, vec.size());

  // evaluate
  const std::vector<Eigen::VectorXd>& outs = split->Evaluate(vec);
  EXPECT_EQ(outs.size(), ind.size());
  for( unsigned int i=0; i<ind.size(); ++i ) {
    EXPECT_NEAR((outs[i]-vec.segment(ind(i), size(i))).norm(), 0.0, 1.0e-10);

    Eigen::MatrixXd jactrue = Eigen::MatrixXd::Zero(size(i), vec.size());
    jactrue.block(0, ind(i), size(i), size(i)) += Eigen::MatrixXd::Identity(size(i), size(i));

    const Eigen::MatrixXd& jac = split->Jacobian(i, 0, vec);
    EXPECT_NEAR((jac-jactrue).norm(), 0.0, 1.0e-10);

    const Eigen::VectorXd sens = Eigen::VectorXd::Random(size(i));
    const Eigen::VectorXd& grad = split->Gradient(i, 0, vec, sens);
    EXPECT_NEAR((grad-jac.transpose()*sens).norm(), 0.0, 1.0e-10);

    const Eigen::VectorXd targ = Eigen::VectorXd::Random(vec.size());
    const Eigen::VectorXd& jacapp = split->ApplyJacobian(i, 0, vec, targ);
    EXPECT_NEAR((jacapp-jac*targ).norm(), 0.0, 1.0e-10);
  }
}

TEST(SplitVectorTests, Extract) {
  // create a vector we want to split
  Eigen::VectorXd vec = Eigen::VectorXd::Random(10);

  // the starting index
  Eigen::VectorXi ind = Eigen::Vector2i(2, 8);

  // the size of each section
  Eigen::VectorXi size = Eigen::Vector2i(4, 2);

  // create the splitter
  auto split = std::make_shared<SplitVector>(ind, size, vec.size());

  // evaluate
  const std::vector<Eigen::VectorXd>& outs = split->Evaluate(vec);
  EXPECT_EQ(outs.size(), ind.size());
  for( unsigned int i=0; i<ind.size(); ++i ) {
    EXPECT_NEAR((outs[i]-vec.segment(ind(i), size(i))).norm(), 0.0, 1.0e-10);

    Eigen::MatrixXd jactrue = Eigen::MatrixXd::Zero(size(i), vec.size());
    jactrue.block(0, ind(i), size(i), size(i)) += Eigen::MatrixXd::Identity(size(i), size(i));

    const Eigen::MatrixXd& jac = split->Jacobian(i, 0, vec);
    EXPECT_NEAR((jac-jactrue).norm(), 0.0, 1.0e-10);

    const Eigen::VectorXd sens = Eigen::VectorXd::Random(size(i));
    const Eigen::VectorXd& grad = split->Gradient(i, 0, vec, sens);
    EXPECT_NEAR((grad-jac.transpose()*sens).norm(), 0.0, 1.0e-10);

    const Eigen::VectorXd targ = Eigen::VectorXd::Random(vec.size());
    const Eigen::VectorXd& jacapp = split->ApplyJacobian(i, 0, vec, targ);
    EXPECT_NEAR((jacapp-jac*targ).norm(), 0.0, 1.0e-10);
  }
}
