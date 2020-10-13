#include <gtest/gtest.h>

#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"

using namespace muq::Approximation;

TEST(Quadrature, ClenshawCurtis) {

  // polynomial order
  int order = 2;

  // legendre table for polynomial order n=5
  std::vector<double> ptsTable = {-1.0000000000000000,
                                  -0.7071067811865476,
                                  0.0000000000000000,
                                  0.7071067811865475,
                                  1.0000000000000000};

  std::vector<double> wtsTable = {0.0666666666666667,
                                  0.5333333333333333,
                                  0.8000000000000000,
                                  0.5333333333333333,
                                  0.0666666666666667};

  // Create quadrature object
  ClenshawCurtisQuadrature quad;

  // Compute pts and Weights
  quad.Compute(order);
  Eigen::VectorXd pts = quad.Points().transpose();
  Eigen::VectorXd wts = quad.Weights();
  EXPECT_EQ(ptsTable.size(), pts.size());
  EXPECT_EQ(wtsTable.size(), wts.size());

  for (int i=0; i<std::min<int>(wts.size(),wtsTable.size()); i++) {
    EXPECT_NEAR(pts(i), ptsTable[i], 1e-9);
    EXPECT_NEAR(wts(i), wtsTable[i], 1e-9);
  }
}
