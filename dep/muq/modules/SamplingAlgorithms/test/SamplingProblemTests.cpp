#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(SamplingProblem, GaussianTarget) {
  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  Eigen::VectorXd mu(2);
  mu << 1,1;
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // it is standard normal (1D) by default

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  Eigen::VectorXd xc = Eigen::VectorXd::Zero(2);
  auto state = std::make_shared<SamplingState>(xc, 1.0);

  EXPECT_DOUBLE_EQ(dist->LogDensity(xc), problem->LogDensity(0, state, AbstractSamplingProblem::SampleType::Accepted));
}
