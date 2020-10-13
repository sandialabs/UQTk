#include <gtest/gtest.h>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"

#include "MUQ/SamplingAlgorithms/MarkovChain.h"


using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

class MarkovChainTest : public::testing::Test {
protected:

    virtual void SetUp() override {

      L.resize(2,2);
      L << 1.0, 0.0,
	1.0, 2.0;

      samps = L * RandomGenerator::GetNormal(2,numWeightedSamps);
      //weights = (maxRepeat*RandomGenerator::GetUniform(numWeightedSamps)).array().ceil();
      weights = Eigen::VectorXd::Ones(numWeightedSamps);
      numSamps = weights.sum();

      for(int i=0; i<numWeightedSamps; ++i)
        collection.Add(std::make_shared<SamplingState>(Eigen::VectorXd(samps.col(i)), weights(i)));
    }

    virtual void TearDown() override {
    }

    const int maxRepeat = 10;
    int numSamps;
    const int numWeightedSamps = 1e5;
    Eigen::MatrixXd L;

    Eigen::MatrixXd samps;
    Eigen::VectorXd weights;
    MarkovChain collection;
};

TEST_F(MarkovChainTest, Mean)
{
  Eigen::VectorXd mu = collection.Mean(0);

  Eigen::VectorXd trueMu = samps*weights / weights.sum();

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 10*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 10*mcStd);

  EXPECT_NEAR(trueMu(0), mu(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu(1), 1e-13);

  Eigen::VectorXd mu2 = collection.Mean();

  mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(0), 10*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(1), 10*mcStd);

  EXPECT_NEAR(trueMu(0), mu2(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu2(1), 1e-13);
}

// Used for timing comparison
TEST_F(MarkovChainTest, SampMean)
{
  Eigen::VectorXd mu = samps * weights / weights.sum();

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 10*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 10*mcStd);
}

TEST_F(MarkovChainTest, ToMatrix)
{
  Eigen::MatrixXd sampMat = collection.AsMatrix();
  unsigned int size = collection.size();
  EXPECT_EQ(numSamps, size);
  EXPECT_EQ(numSamps,sampMat.cols());

  Eigen::VectorXd unweightedMean = sampMat.rowwise().mean();
  Eigen::VectorXd weightedMean = collection.Mean();

  for(int i=0; i<weightedMean.size(); ++i)
    EXPECT_NEAR(weightedMean(i), unweightedMean(i), 1e-13);
}


TEST_F(MarkovChainTest, ToWeights)
{
  Eigen::VectorXd sampWeights = collection.Weights();
  EXPECT_EQ(numSamps, sampWeights.size());
  for( unsigned int i=0; i<numSamps; ++i ) { EXPECT_DOUBLE_EQ(1.0, sampWeights(i)); }
}

TEST_F(MarkovChainTest, ESS)
{
  int totalSteps = collection.size();
  Eigen::VectorXd ess = collection.ESS();

  EXPECT_LE(ess(0), totalSteps);
  EXPECT_LE(ess(1), totalSteps);

  int numWhite = 1e3;
  Eigen::VectorXd whiteNoise = RandomGenerator::GetNormal(numWhite);
  double singleEss = MarkovChain::SingleComponentESS(whiteNoise);

  EXPECT_NEAR(numWhite, singleEss, 250);
}


TEST_F(MarkovChainTest, Variance)
{
  Eigen::VectorXd var = collection.Variance(0);

  Eigen::VectorXd sampMu = samps * weights / weights.sum();
  Eigen::VectorXd sampVar = (samps.colwise() - sampMu).array().pow(2.0).matrix() * weights / weights.sum();

  EXPECT_NEAR(L(0,0)*L(0,0), var(0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var(1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampVar(0), var(0), 1e-13);
  EXPECT_NEAR(sampVar(1), var(1), 1e-13);

  Eigen::VectorXd var2 = collection.Variance();
  EXPECT_NEAR(L(0,0)*L(0,0), var2(0), 20.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var2(1), 150.0/sqrt(double(numSamps)));

}

TEST_F(MarkovChainTest, Covariance)
{
  Eigen::MatrixXd cov = collection.Covariance(0);

  Eigen::VectorXd sampMu  = samps * weights / weights.sum();
  Eigen::MatrixXd sampCov = (samps.colwise()-sampMu) * (weights/weights.sum()).asDiagonal() * (samps.colwise()-sampMu).transpose();

  Eigen::MatrixXd trueCov = L*L.transpose();

  EXPECT_NEAR(trueCov(0,0), cov(0,0), 50.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov(0,1), 100.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov(1,0), 100.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov(1,1), 150.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov(1,1), 1e-13);

  Eigen::MatrixXd cov2 = collection.Covariance();

  EXPECT_NEAR(trueCov(0,0), cov2(0,0), 50.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov2(0,1), 100.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov2(1,0), 100.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov2(1,1), 150.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov2(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov2(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov2(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov2(1,1), 1e-13);
}
