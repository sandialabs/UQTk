#include <gtest/gtest.h>

#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Modeling;

/// A model to whose input/output pairs we can store (3 input dimensions, 2 output dimensions)
class func : public ModPiece {
public:

  inline func() : ModPiece(Eigen::VectorXi::Constant(1, 3), Eigen::VectorXi::Constant(1, 2)) {}

  inline virtual ~func() {}

private:

  inline virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    const Eigen::VectorXd& in = inputs[0];

    Eigen::VectorXd temp = Eigen::VectorXd::Constant(2, std::numeric_limits<double>::quiet_NaN());
    temp << in(0)*in(1), in(2);

    outputs.resize(1);
    outputs[0] = temp;
  }
};

TEST(FlannCache, CreateCache) {
  // create the model whose input/output pairs we want to store
  auto f = std::make_shared<func>();

  // create a cache
  auto cache = std::make_shared<FlannCache>(f);

  // generate some random inputs
  std::vector<Eigen::VectorXd> inputs(10);
  for( auto it=inputs.begin(); it!=inputs.end(); ++it ) { *it = Eigen::Vector3d::Random(); }

  // make sure the size is equal to the number of points that we added
  EXPECT_EQ(cache->Size(), 0);

  // add the inputs to the cache
  cache->Add(inputs);
  
  // make sure the size is equal to the number of points that we added
  EXPECT_EQ(cache->Size(), inputs.size());
  
  // try to add them again
  cache->Add(inputs);
  
  // make sure the size is equal to the number of points that we added with not repeats
  EXPECT_EQ(cache->Size(), inputs.size());

  // make sure they are there
  for( unsigned int i=0; i<inputs.size(); ++i ) { EXPECT_NEAR((cache->at(i)-inputs[i]).norm(), 0.0, 1.0e-10); }
  
  // remove a point from the cache
  cache->Remove(inputs[0]);

  // make sure the size is equal to the number of points that we added minus the one we removed
  EXPECT_EQ(cache->Size(), inputs.size()-1);
  
  // generate some random inputs
  std::vector<Eigen::VectorXd> more_inputs(10);
  for( auto it=more_inputs.begin(); it!=more_inputs.end(); ++it ) { *it = Eigen::Vector3d::Random(); }
  
  // add more inputs to the cache
  cache->Add(more_inputs);
  
  // make sure the size is equal to the number of points that we added minus the one we removed
  EXPECT_EQ(cache->Size(), inputs.size()+more_inputs.size()-1);
  
  // find the 5 nearest neighbors to zero
  std::vector<Eigen::VectorXd> neighbors;
  std::vector<Eigen::VectorXd> result;
  cache->NearestNeighbors((Eigen::Vector3d)Eigen::Vector3d::Zero(), 5, neighbors, result);

  // check the neighbors
  EXPECT_EQ(neighbors.size(), 5);
  for( auto n : neighbors ) { EXPECT_EQ(n.size(), 3); }
  EXPECT_EQ(result.size(), 5);
  for( auto r : result ) { EXPECT_EQ(r.size(), 2); }
}
