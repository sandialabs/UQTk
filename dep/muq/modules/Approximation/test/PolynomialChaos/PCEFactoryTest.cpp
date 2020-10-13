#include <gtest/gtest.h>

#include "MUQ/Approximation/PolynomialChaos/PCEFactory.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;

TEST(PolynomialChaos, PCEFactory_Legendre1D_Gauss)
{
  auto legendre = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(legendre);
  auto quadOrders = std::make_shared<MultiIndex>(1);
  quadOrders->SetValue(0,18);

  PCEFactory factory({quad1d}, quadOrders, {legendre});

  std::vector<Eigen::VectorXd> const& quadPts = factory.QuadPts();

  std::vector<Eigen::VectorXd> quadEvals(quadPts.size());
  for(unsigned int i=0; i<quadPts.size(); ++i){
    quadEvals.at(i).resize(2);
    quadEvals.at(i)(0) = legendre->BasisEvaluate(15,quadPts.at(i)(0));
    quadEvals.at(i)(1) = legendre->BasisEvaluate(11,quadPts.at(i)(0));
  }

  std::shared_ptr<PolynomialChaosExpansion> pce = factory.Compute(quadEvals);
  Eigen::MatrixXd pceCoeffs = pce->GetCoeffs();

  // Make sure the coefficients we expect to be zero are zero
  for(unsigned int i=0; i<pceCoeffs.cols(); ++i){
    if(i!=15){
      EXPECT_NEAR(0.0, pceCoeffs(0,i), 6e-14);
    }else{
      EXPECT_NEAR(1.0, pceCoeffs(0,i), 6e-14);
    }
  }

  for(unsigned int i=0; i<pceCoeffs.cols(); ++i){
    if(i!=11){
      EXPECT_NEAR(0.0, pceCoeffs(1,i), 6e-14);
    }else{
      EXPECT_NEAR(1.0, pceCoeffs(1,i), 6e-14);
    }
  }

  // Make sure the PCE evaluation matches what we'd expect
  Eigen::VectorXd testPt(1);
  testPt(0) = 0.1;

  Eigen::VectorXd pceOut = pce->Evaluate(testPt).at(0);
  Eigen::VectorXd trueOut(2);
  trueOut << legendre->BasisEvaluate(15, testPt(0)), legendre->BasisEvaluate(11, testPt(0));

  EXPECT_NEAR(trueOut(0), pceOut(0),5e-15);
  EXPECT_NEAR(trueOut(1), pceOut(1),5e-15);

}


TEST(PolynomialChaos, PCEFactory_Legendre2D_Gauss)
{
  auto legendre = std::make_shared<Legendre>();
  auto quad1d = std::make_shared<GaussQuadrature>(legendre);
  auto quadOrders = std::make_shared<MultiIndex>(2);
  quadOrders->SetValue(0,10);
  quadOrders->SetValue(1,10);

  PCEFactory factory({quad1d,quad1d}, quadOrders, {legendre,legendre});

  std::vector<Eigen::VectorXd> const& quadPts = factory.QuadPts();

  std::vector<Eigen::VectorXd> quadEvals(quadPts.size());
  for(unsigned int i=0; i<quadPts.size(); ++i){
    quadEvals.at(i).resize(2);
    quadEvals.at(i)(0) = legendre->BasisEvaluate(6,quadPts.at(i)(0))*legendre->BasisEvaluate(7,quadPts.at(i)(1));
    quadEvals.at(i)(1) = legendre->BasisEvaluate(1,quadPts.at(i)(0))*legendre->BasisEvaluate(5,quadPts.at(i)(1));
  }

  std::shared_ptr<PolynomialChaosExpansion> pce = factory.Compute(quadEvals);
  Eigen::MatrixXd pceCoeffs = pce->GetCoeffs();

  // Make sure the PCE evaluation matches what we'd expect
  Eigen::VectorXd testPt(2);
  testPt << 0.1, -0.1;

  Eigen::VectorXd pceOut = pce->Evaluate(testPt).at(0);
  Eigen::VectorXd trueOut(2);
  trueOut << legendre->BasisEvaluate(6, testPt(0)) * legendre->BasisEvaluate(7,testPt(1)),
             legendre->BasisEvaluate(1, testPt(0)) * legendre->BasisEvaluate(5,testPt(1));

  EXPECT_NEAR(trueOut(0), pceOut(0),5e-15);
  EXPECT_NEAR(trueOut(1), pceOut(1),5e-15);

  testPt << 0.8, -0.5;
  pceOut = pce->Evaluate(testPt).at(0);
  trueOut << legendre->BasisEvaluate(6, testPt(0)) * legendre->BasisEvaluate(7,testPt(1)),
             legendre->BasisEvaluate(1, testPt(0)) * legendre->BasisEvaluate(5,testPt(1));

  EXPECT_NEAR(trueOut(0), pceOut(0),5e-15);
  EXPECT_NEAR(trueOut(1), pceOut(1),5e-15);
}
