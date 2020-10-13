#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/UniformBox.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(ImportanceSamplingTests, Setup) {
  // the number of samples
  const unsigned int n = 250000;

  pt::ptree pt;
  pt.put<unsigned int>("ImportanceSampling.NumSamples", n);

  // create a uniform distribution---the sampling problem is built around characterizing this distribution
  auto dist = std::make_shared<UniformBox>(Eigen::RowVector2d(-0.5, 1.0));
  auto target = dist->AsDensity();

  // create a Gaussian distribution---this is the biasing distribution
  auto bias = std::make_shared<Gaussian>(Eigen::VectorXd::Constant(1, 0.25));

  // create an instance of importance sampling
  auto is = std::make_shared<ImportanceSampling>(target, bias, pt.get_child("ImportanceSampling"));

  // generate the samples
  std::shared_ptr<SampleCollection> samps = is->Run();

  // make sure the mean matches
  const Eigen::VectorXd& mean = samps->Mean();
  EXPECT_EQ(mean.size(), 1);
  EXPECT_NEAR(mean(0), 0.25, 1.0e-2);
}
