
#include "gtest/gtest.h"

#include "MUQ/Utilities/RandomGenerator.h"


using namespace muq::Utilities;

TEST(UtilitiesRandomGenerator, UniformReal)
{
  // generate a bunch of uniformly distributed samples and keep track of the mean
  int N       = 100000;
  double mean = 0;

  for (int i = 0; i < N; ++i) {
    mean += RandomGenerator::GetUniform();
  }

  mean /= N;

  EXPECT_NEAR(mean, 0.5, 1e-2);
}

TEST(UtilitiesRandomGenerator, Gaussian)
{
  // generate a bunch of uniformly distributed samples and keep track of the mean
  int N       = 100000;
  double mean = 0;

  for (int i = 0; i < N; ++i) {
    mean += RandomGenerator::GetNormal();
  }

  mean /= N;

  EXPECT_NEAR(mean, 0.0, 2e-2);
}

TEST(UtilitiesRandomGenerator, UniformRealInt)
{
  // generate a bunch of uniformly distributed samples and keep track of the mean
  int N       = 100000;
  double mean = 0;

  for (int i = 0; i < N; ++i) {
    mean += RandomGenerator::GetUniformInt(1, 15);
  }

  mean /= N;

  EXPECT_NEAR(mean, 8, 1e-1);
}

TEST(UtilitiesRandomGenerator, SeedResetTest)
{
  // set the seed
  RandomGenerator::SetSeed(3001);

  // generate 3 normal samples
  double before[5];
  double after[5];

  before[0] = RandomGenerator::GetNormal();
  before[1] = RandomGenerator::GetNormal();
  before[2] = RandomGenerator::GetNormal();
  before[3] = RandomGenerator::GetNormal();
  before[4] = RandomGenerator::GetUniform();

  // reset the seed
  RandomGenerator::SetSeed(3001);
  after[0] = RandomGenerator::GetNormal();
  after[1] = RandomGenerator::GetNormal();
  after[2] = RandomGenerator::GetNormal();
  after[3] = RandomGenerator::GetUniform();
  after[4] = RandomGenerator::GetNormal();

  // compare the before and after
  EXPECT_EQ(before[0], after[0]);
  EXPECT_EQ(before[1], after[1]);
  EXPECT_EQ(before[2], after[2]);
  EXPECT_NE(before[3], after[4]); //check that order matters when we interleave these
  EXPECT_NE(before[4], after[3]);
}

TEST(UtilitiesRandomGenerator, SeedTempSet)
{
	  RandomGenerator::SetSeed(3001);

  // generate 3 normal samples
  double before[5];
    double during[5];
  double after[5];

  before[0] = RandomGenerator::GetNormal();
  before[1] = RandomGenerator::GetNormal();
  before[2] = RandomGenerator::GetNormal();
  before[3] = RandomGenerator::GetNormal();
  before[4] = RandomGenerator::GetUniform();

  {
	  RandomGeneratorTemporarySetSeed(5001);
	    during[0] = RandomGenerator::GetNormal();
  during[1] = RandomGenerator::GetNormal();
  during[2] = RandomGenerator::GetNormal();
  during[3] = RandomGenerator::GetNormal();
  during[4] = RandomGenerator::GetUniform();
  }

  // reset the seed
  RandomGenerator::SetSeed(3001);
  after[0] = RandomGenerator::GetNormal();
  after[1] = RandomGenerator::GetNormal();
  after[2] = RandomGenerator::GetNormal();
  after[3] = RandomGenerator::GetUniform();
  after[4] = RandomGenerator::GetNormal();

  // compare the before and after
  EXPECT_EQ(before[0], after[0]);
  EXPECT_EQ(before[1], after[1]);
  EXPECT_EQ(before[2], after[2]);
  EXPECT_NE(before[3], after[4]); //check that order matters when we interleave these
  EXPECT_NE(before[4], after[3]);

    EXPECT_NE(before[0], during[0]);
  EXPECT_NE(before[1], during[1]);
  EXPECT_NE(before[2], during[2]);
  EXPECT_NE(before[3], during[3]); //check that order matters when we interleave these
  EXPECT_NE(before[4], during[4]);
}

TEST(UtilitiesRandomGenerator, Discrete)
{
  Eigen::VectorXd probs(3);
  probs << 0.2, 0.7, 0.1;

  for(int i=0; i<100; ++i){
    int samp = RandomGenerator::GetDiscrete(probs);
    EXPECT_TRUE((samp >= 0) && (samp < probs.size()));
  }

  const int numSamps = 1000;
  Eigen::VectorXi samps = RandomGenerator::GetDiscrete(probs, numSamps);
  Eigen::VectorXd counts = Eigen::VectorXd::Zero(3);
  for(int i=0; i<numSamps; ++i)
    counts(samps(i)) += 1;

  counts /= double(numSamps);

  EXPECT_NEAR(probs(0), counts(0), 5e-2);
  EXPECT_NEAR(probs(1), counts(1), 5e-2);
  EXPECT_NEAR(probs(2), counts(2), 5e-2);
}

TEST(UtilitiesRandomGenerator, Gamma)
{
  const double alpha = 1.5;
  const double beta = 0.5;
  Eigen::VectorXd samps  = RandomGenerator::GetGamma(alpha, beta, 1e4);

  double mu = samps.mean();
  EXPECT_NEAR(alpha*beta, mu, 1e-2);
}


TEST(UtilitiesRandomGenerator, StoreRng)
{
	auto storedRng = RandomGenerator::CopyGenerator();


  // generate 3 normal samples
  double before[5];
    double during[5];
  double after[5];

  before[0] = RandomGenerator::GetNormal();
  before[1] = RandomGenerator::GetNormal();
  before[2] = RandomGenerator::GetNormal();
  before[3] = RandomGenerator::GetNormal();
  before[4] = RandomGenerator::GetUniform();

  {
	    during[0] = RandomGenerator::GetNormal();
  during[1] = RandomGenerator::GetNormal();
  during[2] = RandomGenerator::GetNormal();
  during[3] = RandomGenerator::GetNormal();
  during[4] = RandomGenerator::GetUniform();
  }

  // reset the seed
  RandomGenerator::SetGenerator(storedRng);
  after[0] = RandomGenerator::GetNormal();
  after[1] = RandomGenerator::GetNormal();
  after[2] = RandomGenerator::GetNormal();
  after[3] = RandomGenerator::GetUniform();
  after[4] = RandomGenerator::GetNormal();

  // compare the before and after
  EXPECT_EQ(before[0], after[0]);
  EXPECT_EQ(before[1], after[1]);
  EXPECT_EQ(before[2], after[2]);
  EXPECT_NE(before[3], after[4]); //check that order matters when we interleave these
  EXPECT_NE(before[4], after[3]);

    EXPECT_NE(before[0], during[0]);
  EXPECT_NE(before[1], during[1]);
  EXPECT_NE(before[2], during[2]);
  EXPECT_NE(before[3], during[3]); //check that order matters when we interleave these
  EXPECT_NE(before[4], during[4]);
}
