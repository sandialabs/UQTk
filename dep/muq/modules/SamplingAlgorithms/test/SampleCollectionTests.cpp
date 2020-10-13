#include <gtest/gtest.h>

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/HDF5/HDF5File.h"

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

class SampleCollectionTest : public::testing::Test {
protected:

    virtual void SetUp() override {

      L.resize(2,2);
      L << 1.0, 0.0,
           1.0, 2.0;

      samps = L * RandomGenerator::GetNormal(2,numSamps);
      weights = RandomGenerator::GetUniform(numSamps);
      weights /= weights.sum();

      for(int i=0; i<numSamps; ++i) {
	auto state = std::make_shared<SamplingState>(Eigen::VectorXd(samps.col(i)), weights(i));
	state->meta["id"] = i;
	state->meta["x norm"] = samps.col(i).norm();
	state->meta["vec2"] = (Eigen::Vector2d)(i*Eigen::Vector2d::Ones());
	state->meta["vec3"] = (Eigen::Vector3d)(i*Eigen::Vector3d::Ones());
	state->meta["vec4"] = (Eigen::Vector4d)(i*Eigen::Vector4d::Ones());
	state->meta["vecX"] = (Eigen::VectorXd)(i*Eigen::VectorXd::Ones(5));

        collection.Add(state);
      }

    }

    virtual void TearDown() override {
    }

    const int numSamps = 1e4;
    Eigen::MatrixXd L;

    Eigen::MatrixXd samps;
    Eigen::VectorXd weights;
    SampleCollection collection;
};

TEST_F(SampleCollectionTest, Mean)
{
  Eigen::VectorXd mu = collection.Mean(0);

  Eigen::VectorXd trueMu = samps*weights;

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 3*mcStd);

  EXPECT_NEAR(trueMu(0), mu(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu(1), 1e-13);

  Eigen::VectorXd mu2 = collection.Mean();

  mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu2(1), 3*mcStd);

  EXPECT_NEAR(trueMu(0), mu2(0), 1e-13);
  EXPECT_NEAR(trueMu(1), mu2(1), 1e-13);
}

// Used for timing comparison
TEST_F(SampleCollectionTest, SampMean)
{
  Eigen::VectorXd mu = samps * weights;

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 3*mcStd);
}

class Func : public ModPiece {
public:
  Func() : ModPiece(Eigen::VectorXi::Constant(1, 2), Eigen::VectorXi::Constant(1, 2)) {}

  virtual ~Func() = default;

private:

  virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    outputs.resize(1);
    outputs[0] = inputs[0];
  }
};

TEST_F(SampleCollectionTest, ExpectedValue)
{
  auto f = std::make_shared<Func>();

  Eigen::VectorXd mu = collection.ExpectedValue(f);

  double mcStd = L(0,0)*L(0,0)/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(0), 3*mcStd);
  mcStd = (L(1,0)*L(1,0) + L(1,1)*L(1,1))/sqrt(double(numSamps));
  EXPECT_NEAR(0.0, mu(1), 3*mcStd);
}

TEST_F(SampleCollectionTest, ToMatrix)
{
  Eigen::MatrixXd sampMat = collection.AsMatrix();

  for(int j=0; j<samps.cols(); ++j){
    for(int i=0; i<samps.rows(); ++i)
      EXPECT_DOUBLE_EQ(samps(i,j), sampMat(i,j));
  }
}

TEST_F(SampleCollectionTest, ToWeights)
{
  Eigen::VectorXd sampWeights = collection.Weights();
  for(int i=0; i<weights.size(); ++i)
    EXPECT_DOUBLE_EQ(weights(i), sampWeights(i));
}

TEST_F(SampleCollectionTest, ESS)
{
  double ess = collection.ESS()(0);
  EXPECT_NEAR(std::pow(weights.sum(),2.0) / weights.array().square().sum(), ess, 5e-11);
}


TEST_F(SampleCollectionTest, Variance)
{
  Eigen::VectorXd var = collection.Variance(0);

  Eigen::VectorXd sampMu = samps * weights;
  Eigen::VectorXd sampVar = (samps.colwise() - sampMu).array().pow(2.0).matrix() * weights;

  EXPECT_NEAR(L(0,0)*L(0,0), var(0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var(1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampVar(0), var(0), 1e-13);
  EXPECT_NEAR(sampVar(1), var(1), 1e-13);

  Eigen::VectorXd var2 = collection.Variance();
  EXPECT_NEAR(L(0,0)*L(0,0), var2(0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(L(1,0)*L(1,0) + L(1,1)*L(1,1), var2(1), 50.0/sqrt(double(numSamps)));

}

TEST_F(SampleCollectionTest, Covariance)
{
  Eigen::MatrixXd cov = collection.Covariance(0);
  std::vector<Eigen::MatrixXd> runCov = collection.RunningCovariance(0);

  Eigen::VectorXd sampMu  = samps * weights;
  Eigen::MatrixXd sampCov = (samps.colwise()-sampMu) * weights.asDiagonal() * (samps.colwise()-sampMu).transpose();

  Eigen::MatrixXd trueCov = L*L.transpose();

  EXPECT_NEAR(trueCov(0,0), runCov[numSamps-1](0,0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), runCov[numSamps-1](0,1), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), runCov[numSamps-1](1,0), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), runCov[numSamps-1](1,1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(trueCov(0,0), cov(0,0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov(0,1), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov(1,0), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov(1,1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov(1,1), 1e-13);

  Eigen::MatrixXd cov2 = collection.Covariance();

  EXPECT_NEAR(trueCov(0,0), cov2(0,0), 5.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(0,1), cov2(0,1), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,0), cov2(1,0), 10.0/sqrt(double(numSamps)));
  EXPECT_NEAR(trueCov(1,1), cov2(1,1), 50.0/sqrt(double(numSamps)));

  EXPECT_NEAR(sampCov(0,0), cov2(0,0), 1e-13);
  EXPECT_NEAR(sampCov(0,1), cov2(0,1), 1e-13);
  EXPECT_NEAR(sampCov(1,0), cov2(1,0), 1e-13);
  EXPECT_NEAR(sampCov(1,1), cov2(1,1), 1e-13);
}

TEST_F(SampleCollectionTest, WriteToFile) {
  const std::string filename = "output.h5";

  // write the collection to file
  collection.WriteToFile(filename);

  auto hdf5file = std::make_shared<HDF5File>(filename);

  const Eigen::MatrixXd samples = hdf5file->ReadMatrix("/samples");
  const Eigen::MatrixXd wghts = hdf5file->ReadMatrix("/weights");
  const Eigen::MatrixXd id = hdf5file->ReadMatrix("/id");
  const Eigen::MatrixXd vec2 = hdf5file->ReadMatrix("/vec2");
  const Eigen::MatrixXd vec3 = hdf5file->ReadMatrix("/vec3");
  const Eigen::MatrixXd vec4 = hdf5file->ReadMatrix("/vec4");
  const Eigen::MatrixXd vecX = hdf5file->ReadMatrix("/vecX");

  for( unsigned int i=0; i<samps.cols(); ++i ) {
    EXPECT_DOUBLE_EQ(i, id(i));
    EXPECT_DOUBLE_EQ((samples.col(i)-samps.col(i)).norm(), 0.0);
    EXPECT_DOUBLE_EQ(weights(i), wghts(i));

    EXPECT_DOUBLE_EQ((vec2.col(i)-i*Eigen::Vector2d::Ones()).norm(), 0.0);
    EXPECT_DOUBLE_EQ((vec3.col(i)-i*Eigen::Vector3d::Ones()).norm(), 0.0);
    EXPECT_DOUBLE_EQ((vec4.col(i)-i*Eigen::Vector4d::Ones()).norm(), 0.0);
    EXPECT_DOUBLE_EQ((vecX.col(i)-i*Eigen::VectorXd::Ones(5)).norm(), 0.0);
  }

  std::remove(filename.c_str());
}
