#include <gtest/gtest.h>

#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;

TEST(PolynomialChaos, PCE_Moments)
{
  const int outDim = 2;
  const int inputDim = 2;
  const int order = 3;

  auto multis = MultiIndexFactory::CreateTotalOrder(inputDim, order);
  auto poly = std::make_shared<Legendre>();

  // Set up a PCE with specific coefficients
  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(outDim, multis->Size());
  Eigen::RowVectorXi multiVec(inputDim);
  multiVec << 1, 2;
  unsigned int ind = multis->MultiToIndex(std::make_shared<MultiIndex>(multiVec));
  coeffs(0,0) = 0.45;
  coeffs(0,ind) = 1.0;
  coeffs(1,ind) = 0.5;

  multiVec << 1, 1;
  ind = multis->MultiToIndex(std::make_shared<MultiIndex>(multiVec));
  coeffs(1,0) = -0.1;
  coeffs(1,ind) = 1.0;

  auto pce = std::make_shared<PolynomialChaosExpansion>(poly, multis, coeffs);

  // Compare the pce mean, variance, and covariance to the analytical values
  Eigen::VectorXd pceMean = pce->Mean();
  EXPECT_DOUBLE_EQ(coeffs(0,0), pceMean(0));
  EXPECT_DOUBLE_EQ(coeffs(1,0), pceMean(1));

  Eigen::MatrixXd trueCov(2,2);
  trueCov << 1.0/15.0, 1.0/30.0,
             1.0/30.0, 23.0/180.0;

  Eigen::MatrixXd pceCov = pce->Covariance();
  EXPECT_DOUBLE_EQ(trueCov(0,0), pceCov(0,0));
  EXPECT_DOUBLE_EQ(trueCov(0,1), pceCov(0,1));
  EXPECT_DOUBLE_EQ(trueCov(1,0), pceCov(1,0));
  EXPECT_DOUBLE_EQ(trueCov(1,1), pceCov(1,1));

  Eigen::VectorXd pceVar = pce->Variance();
  EXPECT_DOUBLE_EQ(trueCov(0,0), pceVar(0));
  EXPECT_DOUBLE_EQ(trueCov(1,1), pceVar(1));
}

TEST(PolynomialChaos, PCE_Sensitivities)
{
  const int outDim = 1;
  const int inputDim = 2;
  const int order = 4;

  auto multis = MultiIndexFactory::CreateTotalOrder(inputDim, order);
  auto poly = std::make_shared<Legendre>();

  // Set up a PCE with specific coefficients
  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(outDim, multis->Size());
  Eigen::RowVectorXi multiVec(inputDim);
  multiVec << 3, 0;
  unsigned int ind = multis->MultiToIndex(std::make_shared<MultiIndex>(multiVec));
  coeffs(0,0) = 0.0;
  coeffs(0,ind) = 1.0;

  multiVec << 1, 1;
  ind = multis->MultiToIndex(std::make_shared<MultiIndex>(multiVec));
  coeffs(0,ind) = 1.0;

  auto pce = std::make_shared<PolynomialChaosExpansion>(poly, multis, coeffs);

  Eigen::VectorXd pceVar = pce->Variance();
  double trueVar = 16.0/63.0;
  EXPECT_DOUBLE_EQ(trueVar, pceVar(0));

  Eigen::MatrixXd pceMain = pce->MainSensitivity();
  Eigen::MatrixXd trueMain(1,2);
  trueMain << (1.0/7.0)/trueVar, 0.0;
  EXPECT_DOUBLE_EQ(trueMain(0),pceMain(0));
  EXPECT_DOUBLE_EQ(trueMain(1),pceMain(1));

  //trueMain << , 0.0;
  Eigen::MatrixXd pceTotal = pce->TotalSensitivity();
  Eigen::MatrixXd trueTotal(1,2);
  trueTotal << 1.0, 1.0-trueMain(0);
  EXPECT_DOUBLE_EQ(trueTotal(0), pceTotal(0));
  EXPECT_DOUBLE_EQ(trueTotal(1), pceTotal(1));
}

TEST(PolynomialChaos, PCE_WeightedSum)
{
  const int outDim = 1;
  const int inputDim = 2;
  const int order = 3;

  auto multis = MultiIndexFactory::CreateTotalOrder(inputDim, order);
  auto poly = std::make_shared<ProbabilistHermite>();

  std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions(3);

  // create an expansion with all coefficients equal to 1
  expansions.at(0) = std::make_shared<PolynomialChaosExpansion>(poly, multis, Eigen::MatrixXd::Ones(outDim,multis->Size()));

  // create an expansion with all coefficients equal to 2
  expansions.at(1) = std::make_shared<PolynomialChaosExpansion>(poly, multis, 2.0*Eigen::MatrixXd::Ones(outDim,multis->Size()));

  // create an expansion with all coefficients equal to 3
  expansions.at(2) = std::make_shared<PolynomialChaosExpansion>(poly, multis, 3.0*Eigen::MatrixXd::Ones(outDim,multis->Size()));

  Eigen::VectorXd weights(3);
  weights << 1.0, 1.0, 0.0;
  auto sumExpansion = PolynomialChaosExpansion::ComputeWeightedSum(expansions,weights);

  Eigen::VectorXd testPt = Eigen::VectorXd::Random(inputDim);

  Eigen::VectorXd trueEval = weights(0)*expansions.at(0)->Evaluate(testPt).at(0) + weights(1)*expansions.at(1)->Evaluate(testPt).at(0);
  EXPECT_NEAR(trueEval(0),sumExpansion->Evaluate(testPt).at(0)(0), 1e-14);

  weights << 0.5, 0.25, 1.0;
  sumExpansion = PolynomialChaosExpansion::ComputeWeightedSum(expansions,weights);
  trueEval = weights(0)*expansions.at(0)->Evaluate(testPt).at(0) + weights(1)*expansions.at(1)->Evaluate(testPt).at(0) + weights(2)*expansions.at(2)->Evaluate(testPt).at(0);
  EXPECT_NEAR(trueEval(0),sumExpansion->Evaluate(testPt).at(0)(0), 1e-14);
}
