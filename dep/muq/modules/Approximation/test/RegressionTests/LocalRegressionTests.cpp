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

class LocalRegressionTest : public::testing::Test {
public:

  inline LocalRegressionTest() {
    // the function too approximate
    fn = std::make_shared<func>();

    // set the regressor options
    pt::ptree pt;
    pt.put<unsigned int>("LocalRegression.NumNeighbors", 11);
    pt.put<unsigned int>("LocalRegression.Order", 2);

    // create a local regressor
    reg = std::make_shared<LocalRegression>(fn, pt.get_child("LocalRegression"));

    // generate some random inputs
    std::vector<Eigen::VectorXd> inputs(M);
    for( auto it=inputs.begin(); it!=inputs.end(); ++it ) { *it = Eigen::Vector3d::Random(); }

    // add the random input points to the cache
    reg->Add(inputs);
  }

  inline virtual ~LocalRegressionTest() {}

protected:

  /// The function to approximate
  std::shared_ptr<func> fn;

  /// The local regressor
  std::shared_ptr<LocalRegression> reg;

  // The number of points we added
  const unsigned int M = 25;
};

TEST_F(LocalRegressionTest, Basic) {
  // check the size
  EXPECT_EQ(reg->CacheSize(), M);

  // the input point
  const Eigen::VectorXd input = Eigen::Vector3d::Random();

  // evaluate the local polynomial approximation
  const std::vector<Eigen::VectorXd>& result = reg->Evaluate(input);

  // evaluate the truth
  const std::vector<Eigen::VectorXd>& truth = fn->Evaluate(input);

  // the regression and the truth are the same---approximating a quadratic with a quardratic
  EXPECT_NEAR((truth[0]-result[0]).norm(), 0.0, 1.0e-10);
}

TEST_F(LocalRegressionTest, Poisedness) {
  // check the size
  EXPECT_EQ(reg->CacheSize(), M);

  for( unsigned int i=0; i<10; ++i ) {
    // the input point
    const Eigen::VectorXd input = Eigen::Vector3d::Random();

    std::vector<Eigen::VectorXd> neighbors;
    std::vector<Eigen::VectorXd> results;
    reg->NearestNeighbors(input, neighbors, results);
    EXPECT_TRUE(neighbors.size()==results.size());

    // get the poisedness constant
    std::tuple<Eigen::VectorXd, double, unsigned int> lambda = reg->PoisednessConstant(input, neighbors);
    EXPECT_TRUE(std::get<2>(lambda)<neighbors.size());

    const Eigen::VectorXd& newResult = reg->Add(std::get<0>(lambda));
    EXPECT_EQ(reg->CacheSize(), M+i+1);
  }
}
