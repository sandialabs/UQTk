#include "MUQ/config.h"

#if MUQ_HAS_MPI

#include <gtest/gtest.h>

#include <Eigen/Core>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/SamplingAlgorithms/DistributedCollection.h"

using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

class DistributedCollectionTest : public::testing::Test {
protected:

  inline virtual void SetUp() override {
    L.resize(2,2);
    L << 1.0, 0.0,
      1.0, 2.0;

    // the default communicator is MPI_COMM_WORLD
    auto comm = std::make_shared<parcer::Communicator>();
    rank = comm->GetRank();
    nproc = comm->GetSize();

    // create the local collection
    auto local = std::make_shared<SampleCollection>();

    for( unsigned int i=0; i<numSamps; ++i ) {

      auto state = std::make_shared<SamplingState>((Eigen::VectorXd)(L*RandomGenerator::GetNormal(2)), (double)RandomGenerator::GetUniform());
      state->meta["rank"] = rank;
      state->meta["id"] = i;

      local->Add(state);
    }

    // create the global collection
    collection = std::make_shared<DistributedCollection>(local, comm);

  }

  inline virtual void TearDown() override {

    int cnt = 0;
    for( unsigned int i=0; i<numSamps; ++i ) {

      auto state = collection->LocalAt(i);
      EXPECT_EQ(boost::any_cast<unsigned int const>(state->meta["id"]), i);
      EXPECT_EQ(boost::any_cast<int const>(state->meta["rank"]), rank);
    }

    auto state0 = collection->GlobalAt(0);
    EXPECT_EQ(boost::any_cast<int const>(state0->meta["rank"]), 0);
    EXPECT_EQ(boost::any_cast<unsigned int const>(state0->meta["id"]),0);

    state0 = collection->GlobalAt(numSamps);
    EXPECT_EQ(boost::any_cast<int const>(state0->meta["rank"]), 1);
    EXPECT_EQ(boost::any_cast<unsigned int const>(state0->meta["id"]),0);
  }

  const int numSamps = 5e5;
  Eigen::MatrixXd L;

  int rank;
  int nproc;

  std::shared_ptr<DistributedCollection> collection;
};

TEST_F(DistributedCollectionTest, SizeTest) {
  EXPECT_EQ(collection->LocalSize(), numSamps);
  EXPECT_EQ(collection->GlobalSize(), nproc*numSamps);
  EXPECT_EQ(collection->size(), 2*numSamps);
}

TEST_F(DistributedCollectionTest, CentralMoment) {
  EXPECT_NEAR(collection->LocalCentralMoment(1).norm(), 0.0, 1.0e-2);
  EXPECT_NEAR(collection->GlobalCentralMoment(1).norm(), 0.0, 1.0e-2);
  EXPECT_NEAR(collection->CentralMoment(1).norm(), 0.0, 1.0e-2);
}

TEST_F(DistributedCollectionTest, Mean) {
  EXPECT_NEAR(collection->LocalMean().norm(), 0.0, 1.0e-1);
  EXPECT_NEAR(collection->GlobalMean().norm(), 0.0, 1.0e-1);
  EXPECT_NEAR(collection->Mean().norm(), 0.0, 1.0e-1);
}

TEST_F(DistributedCollectionTest, Variance) {
  const Eigen::VectorXd& truth = (L*L.transpose()).diagonal();

  EXPECT_NEAR((collection->LocalVariance()-truth).norm(), 0.0, 1.0e-1);
  EXPECT_NEAR((collection->GlobalVariance()-truth).norm(), 0.0, 1.0e-1);
  EXPECT_NEAR((collection->Variance()-truth).norm(), 0.0, 1.0e-1);
}

TEST_F(DistributedCollectionTest, Covariance) {
  const Eigen::MatrixXd& truth = L*L.transpose();

  EXPECT_NEAR((collection->LocalCovariance()-truth).norm(), 0.0, 1.0e-1);
  EXPECT_NEAR((collection->GlobalCovariance()-truth).norm(), 0.0, 1.0e-1);
  EXPECT_NEAR((collection->Covariance()-truth).norm(), 0.0, 1.0e-1);
}

TEST_F(DistributedCollectionTest, ESS) {
  const Eigen::VectorXd& less = collection->LocalESS();
  const Eigen::VectorXd& gess = collection->GlobalESS();
  const Eigen::VectorXd& ess = collection->ESS();

  EXPECT_EQ(less.size(), ess.size());
  EXPECT_EQ(gess.size(), ess.size());

  for( unsigned int i=0; i<ess.size(); ++i ) {
    EXPECT_TRUE(less(i)>0);
    EXPECT_TRUE(gess(i)>0);
    EXPECT_TRUE(ess(i)>0);
  }
}

TEST_F(DistributedCollectionTest, AsMatrix) {
  const Eigen::MatrixXd& localMat = collection->AsLocalMatrix();
  EXPECT_EQ(localMat.rows(), 2);
  EXPECT_EQ(localMat.cols(), numSamps);

  const Eigen::MatrixXd& globalMat = collection->AsGlobalMatrix();
  EXPECT_EQ(globalMat.rows(), 2);
  EXPECT_EQ(globalMat.cols(), nproc*numSamps);

  const Eigen::MatrixXd& mat = collection->AsMatrix();
  EXPECT_EQ(mat.rows(), 2);
  EXPECT_EQ(mat.cols(), nproc*numSamps);
}

TEST_F(DistributedCollectionTest, Weights) {
  const Eigen::VectorXd& localWeights = collection->LocalWeights();
  EXPECT_EQ(localWeights.size(), numSamps);

  const Eigen::MatrixXd& globalWeights = collection->GlobalWeights();
  EXPECT_EQ(globalWeights.size(), nproc*numSamps);

  const Eigen::MatrixXd& weights = collection->Weights();
  EXPECT_EQ(weights.size(), nproc*numSamps);
}

#endif // end MUQ_HAS_MPI
