#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/ParallelAMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

TEST(MCMC, MHKernel_MHProposal) {
  const unsigned int N = 1e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "MHProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 0.5); // the variance of the isotropic MH proposal

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = mu;

  auto comm = std::make_shared<parcer::Communicator>();

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem, comm);

  std::shared_ptr<SampleCollection> localSamps = mcmc->Run(start);
  auto samps = std::make_shared<DistributedCollection>(localSamps, comm);
  EXPECT_EQ(comm->GetSize()*pt.get<int>("MyMCMC.NumSamples"), samps->size());
  
  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);
  
  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}

TEST(MCMC, MHKernel_ParallelAMProposal) {
  const unsigned int N = 1.0e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("MyMCMC.NumSamples", N); // number of Monte Carlo samples
  pt.put("MyMCMC.PrintLevel",0);
  pt.put("MyMCMC.KernelList", "Kernel1"); // the transition kernel
  pt.put("MyMCMC.Kernel1.Method","MHKernel");
  pt.put("MyMCMC.Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("MyMCMC.Kernel1.MyProposal.Method", "ParallelAMProposal");
  pt.put("MyMCMC.Kernel1.MyProposal.ProposalVariance", 1.0); // the variance of the isotropic MH proposal
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptSteps",200);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptStart",2000);
  pt.put("MyMCMC.Kernel1.MyProposal.AdaptScale",1.0);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity(); // standard normal Gaussian

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  // starting point
  const Eigen::VectorXd start = Eigen::VectorXd::Random(2);

  auto comm = std::make_shared<parcer::Communicator>();

  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt.get_child("MyMCMC"), problem, comm);

  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<ParallelAMProposal> proposalAM = std::dynamic_pointer_cast<ParallelAMProposal>(proposalBase);
  EXPECT_TRUE(proposalAM);

  std::shared_ptr<SampleCollection> localSamps = mcmc->Run(start);
  auto samps = std::make_shared<DistributedCollection>(localSamps, comm);
  EXPECT_EQ(comm->GetSize()*pt.get<int>("MyMCMC.NumSamples"), samps->size());

  Eigen::VectorXd mean = samps->Mean();
  EXPECT_NEAR(mu(0), mean(0), 1e-1);
  EXPECT_NEAR(mu(1), mean(1), 1e-1);

  Eigen::MatrixXd cov = samps->Covariance();
  EXPECT_NEAR(1.0, cov(0,0), 1e-1);
  EXPECT_NEAR(0.0, cov(0,1), 1e-1);
  EXPECT_NEAR(0.0, cov(1,0), 1e-1);
  EXPECT_NEAR(1.0, cov(1,1), 1e-1);
}
