#include <gtest/gtest.h>

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

TEST(MCMC, CrankNicolson) {

  const unsigned int N = 5e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("BurnIn", 1e4);
  pt.put("PrintLevel",0);
  pt.put("KernelList", "Kernel1"); // the transition kernel
  pt.put("Kernel1.Method","MHKernel");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "CrankNicolsonProposal");
  pt.put("Kernel1.MyProposal.Beta", 0.8);
  pt.put("Kernel1.MyProposal.PriorNode", "Prior"); // The node in the WorkGraph containing the prior density

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = 1.5*Eigen::VectorXd::Ones(2);
  const Eigen::MatrixXd priorCov = Eigen::MatrixXd::Identity(2,2);
  auto priorDist = std::make_shared<Gaussian>(mu, priorCov);

  Eigen::VectorXd obs(2);
  obs << 1.0, 2.0;
  Eigen::MatrixXd obsCov(2,2);
  obsCov << 1.0, 0.05,
            0.05, 1.0;
  auto obsDist = std::make_shared<Gaussian>(obs, obsCov);

  auto graph = std::make_shared<WorkGraph>();
  graph->AddNode(std::make_shared<IdentityOperator>(2), "Parameters");
  graph->AddNode(priorDist->AsDensity(), "Prior");
  graph->AddNode(obsDist->AsDensity(), "Likelihood");
  graph->AddEdge("Parameters", 0, "Prior", 0);
  graph->AddEdge("Parameters", 0, "Likelihood", 0);

  graph->AddNode(std::make_shared<DensityProduct>(2), "Posterior");
  graph->AddEdge("Prior",0,"Posterior",0);
  graph->AddEdge("Likelihood",0,"Posterior",1);

  auto dens = graph->CreateModPiece("Posterior");

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dens);


  // starting point
  std::vector<Eigen::VectorXd> start(1);
  start.at(0) = mu;

  // Compute the true posterior mean and covariance
  Eigen::MatrixXd predCov = obsCov + priorCov;
  Eigen::VectorXd postMean = mu + priorCov * predCov.selfadjointView<Eigen::Lower>().llt().solve(obs- mu);
  Eigen::MatrixXd postCov = priorCov - priorCov * predCov.selfadjointView<Eigen::Lower>().llt().solve(priorCov);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);

  // Make sure the kernel and proposal are correct
  EXPECT_EQ(1,mcmc->Kernels().size());

  std::shared_ptr<TransitionKernel> kernelBase = mcmc->Kernels().at(0);
  std::shared_ptr<MHKernel> kernelMH = std::dynamic_pointer_cast<MHKernel>(kernelBase);
  EXPECT_TRUE(kernelMH);

  std::shared_ptr<MCMCProposal> proposalBase = kernelMH->Proposal();
  std::shared_ptr<CrankNicolsonProposal> proposalCN = std::dynamic_pointer_cast<CrankNicolsonProposal>(proposalBase);
  EXPECT_TRUE(proposalCN);

  std::shared_ptr<SampleCollection> samps = mcmc->Run(start);

  Eigen::VectorXd sampMean = samps->Mean();
  Eigen::MatrixXd sampCov = samps->Covariance();

  EXPECT_NEAR(postMean(0), sampMean(0), 1e-2);
  EXPECT_NEAR(postMean(1), sampMean(1), 1e-2);

}


TEST(MCMC, CrankNicolsonInGibbs) {

  const unsigned int N = 1e5;

  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("BurnIn", 1e4);
  pt.put("PrintLevel",0);

  pt.put("KernelList", "Kernel1,Kernel2"); // the transition kernel

  pt.put("Kernel1.Method","MHKernel");
  pt.put("Kernel1.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel1.MyProposal.Method", "CrankNicolsonProposal");
  pt.put("Kernel1.MyProposal.Beta", 0.8);
  pt.put("Kernel1.MyProposal.PriorNode", "Prior"); // The node in the WorkGraph containing the prior density

  pt.put("Kernel2.Method","MHKernel");
  pt.put("Kernel2.Proposal", "MyProposal"); // the proposal
  pt.put("Kernel2.MyProposal.Method", "MHProposal");
  pt.put("Kernel2.MyProposal.ProposalVariance", 0.1);

  // create a Gaussian distribution---the sampling problem is built around characterizing this distribution
  const Eigen::VectorXd mu = 1.5*Eigen::VectorXd::Ones(2);
  auto priorDist = std::make_shared<Gaussian>(mu, Gaussian::DiagCovariance);

  Eigen::VectorXd obs(2);
  obs << 1.0, 2.0;

  Eigen::MatrixXd obsCov(2,2);
  obsCov << 1.0, 0.05,
            0.05, 1.0;

  auto obsDist = std::make_shared<Gaussian>(obs, obsCov);

  auto graph = std::make_shared<WorkGraph>();
  graph->AddNode(std::make_shared<IdentityOperator>(2), "Parameters");
  graph->AddNode(priorDist->AsDensity(), "Prior");
  graph->AddNode(obsDist->AsDensity(), "Likelihood");
  graph->AddEdge("Parameters", 0, "Prior", 0);
  graph->AddEdge("Parameters", 0, "Likelihood", 0);

  graph->AddNode(std::make_shared<IdentityOperator>(2), "Test");
  graph->AddEdge("Test", 0, "Prior",1);

  graph->AddNode(std::make_shared<DensityProduct>(2), "Posterior");
  graph->AddEdge("Prior",0,"Posterior",0);
  graph->AddEdge("Likelihood",0,"Posterior",1);

  auto dens = graph->CreateModPiece("Posterior");

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dens);


  // starting point
  std::vector<Eigen::VectorXd> start(2);
  start.at(0) = mu;
  start.at(1) = Eigen::VectorXd::Ones(2);

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<SingleChainMCMC>(pt,problem);

  // Make sure the kernel and proposal are correct
  EXPECT_EQ(2,mcmc->Kernels().size());

  mcmc->Run(start);
}
