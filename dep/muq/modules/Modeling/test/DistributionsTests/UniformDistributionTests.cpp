#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/UniformBox.h"

using namespace muq::Modeling;

class UniformDistributionTests : public::testing::Test {
public:

  inline UniformDistributionTests() {
    const std::pair<double, double> first_bounds(0.0, 1.0);
    const std::pair<double, double> second_bounds(-2.0, 3.0);

    // create a uniform distribution
    uniform1D = std::make_shared<UniformBox>(0.0, 1.0);//first_bounds);

    uniform = std::make_shared<UniformBox>(first_bounds, second_bounds);
    uniform2 = std::make_shared<UniformBox>(0.0, 1.0, -2.0, 3.0);
  }

  inline virtual ~UniformDistributionTests() {}

  std::shared_ptr<UniformBox> uniform, uniform2, uniform1D;

private:
};

TEST_F(UniformDistributionTests, EvaluateLogDensity) {
  // a point inside the domain
  Eigen::VectorXd x0(2);
  x0 << 0.25, 0.0;

  // a point outside the domain
  Eigen::VectorXd x1(2);
  x1 << 0.25, -4.0;

  // a point outside the domain
  Eigen::VectorXd x2(2);
  x2 << 2.8, -4.0;

  // evalute the log-denstion
  double logdens = uniform->LogDensity(x0);
  EXPECT_DOUBLE_EQ(logdens, -log(5.0));

  logdens = uniform2->LogDensity(x0);
  EXPECT_DOUBLE_EQ(logdens, -log(5.0));

  // evalute the log-denstion
  logdens = uniform->LogDensity(x1);
  EXPECT_DOUBLE_EQ(logdens, -std::numeric_limits<double>::infinity());

  // evalute the log-denstion
  logdens = uniform->LogDensity(x2);
  EXPECT_DOUBLE_EQ(logdens, -std::numeric_limits<double>::infinity());
}

TEST_F(UniformDistributionTests, Sample1D) {
  // make sure the expected value is correct
  const unsigned int N = 100;

  for( unsigned int i=0; i<N; ++i ) {

      const double test = uniform1D->Sample()(0);
      EXPECT_TRUE(test>=0.0);
      EXPECT_TRUE(test<=1.0);
  }
}

TEST_F(UniformDistributionTests, Sample) {
  const Eigen::VectorXd samp = uniform->Sample();

  EXPECT_TRUE(samp(0)>=0.0);
  EXPECT_TRUE(samp(0)<=1.0);
  EXPECT_TRUE(samp(1)>=-2.0);
  EXPECT_TRUE(samp(1)<=3.0);

}
