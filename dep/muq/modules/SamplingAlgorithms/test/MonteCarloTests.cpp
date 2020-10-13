#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/MonteCarlo.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(MonteCarlo, GaussianTarget) {
  // create an instance of Monte Carlo
  auto mc = std::make_shared<MonteCarlo>();

  // the number of samples
  const unsigned int N = 1.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put<unsigned int>("SamplingAlgorithm.NumSamples", N); // number of Monte Carlo samples

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<Gaussian>(); // it is standard normal (1D) by default

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // evaluate
  const std::vector<boost::any>& result = mc->Evaluate(pt, problem);
  const std::vector<std::shared_ptr<SamplingState> >& samples = boost::any_cast<std::vector<std::shared_ptr<SamplingState> > const&>(result[0]);
  EXPECT_EQ(samples.size(), N);

  // estimate the mean
  const boost::any mean = mc->FirstMoment();

  EXPECT_NEAR(boost::any_cast<double const>(mean), 0.0, 1.0e-2);
}

