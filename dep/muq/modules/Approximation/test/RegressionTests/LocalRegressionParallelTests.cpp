#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/LocalRegression.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

/// A model to approximate (3 input dimensions, 2 output dimensions)
class func : public ModPiece {
public:

  inline func() : ModPiece(Eigen::VectorXi::Constant(1,3), Eigen::VectorXi::Constant(1,2)) {}

  inline virtual ~func() {}

private:

  inline virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& ins) override {
    outputs.resize(1);
    outputs[0] = (Eigen::VectorXd)Eigen::Vector2d(ins[0](0)*ins[0](1), ins[0](2));
  }
};

TEST(LocalRegressionTest, SharedCache) {
  // the function too approximate
  auto fn = std::make_shared<func>();

  // the communicator
  auto comm = std::make_shared<parcer::Communicator>();
  comm->Barrier();

  // set the regressor options
  pt::ptree pt;
  pt.put<unsigned int>("LocalRegression.NumNeighbors", 11);
  pt.put<unsigned int>("LocalRegression.Order", 2);

  // create a local regressor
  auto reg = std::make_shared<LocalRegression>(fn, pt.get_child("LocalRegression"), comm);

  // add points one at a time
  const unsigned int N = 25;
  for( unsigned int i=0; i<N; ++i ) { reg->Add(Eigen::Vector3d::Random()); }

  // make sure they all have evalauted the function
  comm->Barrier();

  // get any residual messages
  reg->Probe();

  // check the size
  EXPECT_EQ(reg->CacheSize(), comm->GetSize()*N);

  // make sure they all have evalauted the function
  comm->Barrier();

  // add points as a group
  const unsigned int M = 25;
  std::vector<Eigen::VectorXd> inputs(M);
  for( auto it=inputs.begin(); it!=inputs.end(); ++it ) { *it = Eigen::Vector3d::Random(); }

  // add the random input points to the cache
  reg->Add(inputs);

  // make sure they all have evalauted the function
  comm->Barrier();

  // get any residual messages
  reg->Probe();

  // check the size
  EXPECT_EQ(reg->CacheSize(), comm->GetSize()*(N+M));
}
