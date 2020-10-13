#include <gtest/gtest.h>

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"

#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Laguerre.h"
#include "MUQ/Approximation/Polynomials/Jacobi.h"

#include <Eigen/Core>

using namespace muq::Approximation;

TEST(Quadrature, PhysicistHermite) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {-2.020182870456085632929,
                                       -0.9585724646138185071128,
                                       0.0,
                                       0.9585724646138185071128,
                                       2.020182870456085632929};

  std::vector<double> gaussWtsTable = {0.01995324205904591320774,
                                       0.3936193231522411598285,
                                       0.9453087204829418812257,
                                       0.3936193231522411598285,
                                       0.01995324205904591320774};

  // create a PhysicistHermite object
  auto poly = std::make_shared<PhysicistHermite>();

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, ProbabilistHermite) {

  // polynomial order
  int order = 6;

  // legendre table for polynomial order n=7
  std::vector<double> gaussPtsTable = {-3.750439717725742,
                                       -2.366759410734541,
                                       -1.154405394739968,
                                       0.0,
                                       1.154405394739968,
                                       2.366759410734542,
                                       3.750439717725742};

  std::vector<double> gaussWtsTable = {0.1374306216479550e-02,
                                       0.7709667658348317e-01,
                                       0.6018995488855944,
                                       1.145887211259887,
                                       0.6018995488855950,
                                       0.7709667658348326e-01,
                                       0.1374306216479553e-02};

  // create a ProbabilistHermite object
  auto poly = std::make_shared<ProbabilistHermite>();

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, LegendreFromString) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {-0.9061798459386639927976,
                                       -0.5384693101056830910363,
                                       0.0,
                                       0.5384693101056830910363,
                                       0.9061798459386639927976};

  std::vector<double> gaussWtsTable = {0.2369268850561890875143,
                                       0.4786286704993664680413,
                                       0.5688888888888888888889,
                                       0.4786286704993664680413,
                                       0.2369268850561890875143};

  // create a Legendre object
  auto poly = OrthogonalPolynomial::Construct("Legendre");

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}


TEST(Quadrature, Legendre) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {-0.9061798459386639927976,
                                       -0.5384693101056830910363,
                                       0.0,
                                       0.5384693101056830910363,
                                       0.9061798459386639927976};

  std::vector<double> gaussWtsTable = {0.2369268850561890875143,
                                       0.4786286704993664680413,
                                       0.5688888888888888888889,
                                       0.4786286704993664680413,
                                       0.2369268850561890875143};

  // create a Legendre object
  auto poly = std::make_shared<Legendre>();

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, LaguerreDefault) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {0.2635603197181409102031,
                                       1.413403059106516792218,
                                       3.596425771040722081223,
                                       7.085810005858837556922,
                                       12.64080084427578265943};

  std::vector<double> gaussWtsTable = {0.521755610582808652476,
                                       0.398666811083175927454,
                                       0.0759424496817075953877,
                                       0.00361175867992204845446,
                                       2.33699723857762278911e-5};

  // create a Laguerre object
  auto poly = std::make_shared<Laguerre>();

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, Laguerre) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {1.490554945186828158609,
                                       3.581333812903254000336,
                                       6.62699629682360783548,
                                       10.94441800346766931403,
                                       17.35669694161864069155};

  std::vector<double> gaussWtsTable = {1.250983613344565671655,
                                       3.238557188892596223057,
                                       1.390185244072203664795,
                                       0.1190411728383792356206,
                                       0.001232780852255204873499};

  // Laguerre powers
  double alpha = 3.0;

  // create a Laguerre object
  auto poly = std::make_shared<Laguerre>(alpha);

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, JacobiDefault) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {-0.830223896278566929872,
                                       -0.4688487934707142138038,
                                       0.0,
                                       0.4688487934707142138038,
                                       0.830223896278566929872};

  std::vector<double> gaussWtsTable = {0.0860176821228074532418,
                                       0.336839460734335403901,
                                       0.487619047619047619048,
                                       0.336839460734335403901,
                                       0.08601768212280745324181};

  // create a Jacobi object
  auto poly = std::make_shared<Jacobi>();

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}

TEST(Quadrature, Jacobi) {

  // polynomial order
  int order = 4;

  // legendre table for polynomial order n=5
  std::vector<double> gaussPtsTable = {-0.6904577501267610633887,
                                       -0.32651993134900065265,
                                       0.0,
                                       0.4751788706128316398344,
                                       0.792794294644228504216};

  std::vector<double> gaussWtsTable = {0.0274101780663370994418,
                                       0.2129178606036482785255,
                                       0.4390843794439508053922,
                                       0.3222065654722182153781,
                                       0.0650476830805122679291};

  // Jacobi powers
  double alpha = 3.0;
  double beta = 3.0;

  // create a Jacobi object
  auto poly = std::make_shared<Jacobi>(alpha, beta);

  // Create quadrature object
  GaussQuadrature gq(poly);

  // Compute pts and Weights
  gq.Compute(order);
  Eigen::VectorXd gaussPts = gq.Points().transpose();
  Eigen::VectorXd gaussWts = gq.Weights();

  std::cout << gaussPts << std::endl;
  for (int i=0; i<order+1; i++) {

    EXPECT_NEAR(gaussPts(i), gaussPtsTable[i], 1e-9);
    EXPECT_NEAR(gaussWts(i), gaussWtsTable[i], 1e-9);

  }

}
