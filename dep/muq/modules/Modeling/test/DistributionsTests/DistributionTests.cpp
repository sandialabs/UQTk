#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

class ExampleDensity : public Distribution {
public:

  inline ExampleDensity() : Distribution(1) {}

  inline virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
    // get the point where we are evaluating the log density
    const double x = inputs.at(0)(0);

    return -x*x-std::sin(2.0*M_PI*x);
  }

private:
};

class ExampleRV : public Distribution {
public:
    inline ExampleRV() : Distribution(1) {}

    inline virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override {
        double outputVal = 0.1;
        return outputVal*Eigen::VectorXd::Ones(1);
    }

private:
};

TEST(Distribution, EvaluateDensity) {
  // create a distribution that only has a (log) density
  auto dens = std::make_shared<ExampleDensity>();

  // evaluate the density at a point
  Eigen::VectorXd x(1);
  x << 2.5;

  // evaluate the log density
  double logdens = dens->LogDensity(x);

  // make sure we get the density we implemented
  EXPECT_DOUBLE_EQ(logdens, -x(0)*x(0)-std::sin(2.0*M_PI*x(0)));

  std::shared_ptr<Density> densPiece = dens->AsDensity();
  ASSERT_TRUE(densPiece);

  std::vector<Eigen::VectorXd> const& result = densPiece->Evaluate(x);
  double logDens = result.at(0)(0);
  EXPECT_DOUBLE_EQ(logDens, -x(0)*x(0)-std::sin(2.0*M_PI*x(0)));

  // Make sure we can't sample
  EXPECT_THROW(dens->Sample(x), muq::NotImplementedError);
}

TEST(Distribution, EvaluateSample) {

    auto rv = std::make_shared<ExampleRV>();
    Eigen::VectorXd x(1);
    x << 2.5;

    // draw a sample
    Eigen::VectorXd samp = rv->Sample();
    EXPECT_DOUBLE_EQ(0.1, samp(0));

    samp = rv->AsVariable()->Evaluate().at(0);
    EXPECT_DOUBLE_EQ(0.1, samp(0));

    // Make sure we can't evaluate the density
    EXPECT_THROW(rv->LogDensity(x), muq::NotImplementedError);
}
