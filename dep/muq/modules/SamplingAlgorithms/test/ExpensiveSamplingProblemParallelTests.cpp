#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(ExpensiveSamplingProblemTests, GaussianTarget) {
  const unsigned int N = 1000;

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  pt.put("MySamplingProblem.RegressionOptions", "MyRegression");
  pt.put("MySamplingProblem.MyRegression.NumNeighbors", 12);
  pt.put("MySamplingProblem.MyRegression.Order", 3); // approximating the quardatic log-Gaussian with a locally linear function

  pt.put("MySamplingProblem.BetaScale", 1.0);
  pt.put("MySamplingProblem.BetaExponent", 0.9);

  pt.put("MySamplingProblem.FirstLevelLength", 1.0);

  pt.put("MySamplingProblem.GammaScale", 1.0);
  pt.put("MySamplingProblem.GammaExponent", 0.5);
  pt.put("MaximumGammaRefine", 15);
  pt.put("TargetMax", std::exp(dist->LogDensity(mu)));

  pt.put("MySamplingProblem.EtaScale", 2.0);
  pt.put("MySamplingProblem.EtaExponent", 2.0);

  pt.put("MySamplingProblem.LambdaScale", 25.0);

  // the communicator
  auto comm = std::make_shared<parcer::Communicator>();

  // create a sampling problem
  auto problem = std::make_shared<ExpensiveSamplingProblem>(dist, pt.get_child("MySamplingProblem"), comm);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem);

  // run MCMC
  std::shared_ptr<SampleCollection> localSamps = mcmc->Run(start);

  // make sure the number of evaluations is less than the number of steps
  EXPECT_TRUE(problem->CacheSize()>pt.get<unsigned int>("MySamplingProblem.MyRegression.NumNeighbors"));
  EXPECT_TRUE(problem->CacheSize()<N);
}
