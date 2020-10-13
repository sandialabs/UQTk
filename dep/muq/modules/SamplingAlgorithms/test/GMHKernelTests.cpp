#include <gtest/gtest.h>

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/GMHKernel.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

TEST(GMHKernelTest, ProposalTest) {
  // create the density we wish to sample
  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(2);
  auto dist = std::make_shared<Gaussian>(mu)->AsDensity();

  // create a sampling problem
  auto problem = std::make_shared<SamplingProblem>(dist);

  const unsigned int N = 10;
  const unsigned int M = 10;

  // use an adpative Metropolis proposal
  pt::ptree pt;
  pt.put("MyKernel.MyProposal.ProposalVariance", 1.0);
  pt.put("MyKernel.MyProposal.AdaptSteps", 200);
  pt.put("MyKernel.MyProposal.AdaptStart", 2000);
  pt.put("MyKernel.MyProposal.AdaptScale", 1.0);
  pt.put("MyKernel.MyProposal.AdaptScale", 1.0);
  pt.put("MyKernel.MyProposal.Method", "AMProposal");
  pt.put("MyKernel.Proposal", "MyProposal");
  pt.put("MyKernel.NumProposals", N);
  pt.put("MyKernel.NumAccepted", M); // optional: defaults to N

  // create the kernel
  auto kern = std::make_shared<GMHKernel>(pt.get_child("MyKernel"), problem);

  // some dummy state
  auto state = std::make_shared<SamplingState>(Eigen::VectorXd::Random(2), 1.0);

  // propose a bunch of points in the pre step
  kern->PreStep(0, state);
  const Eigen::VectorXd& stationaryAcceptance = kern->StationaryAcceptance();
  EXPECT_EQ(stationaryAcceptance.size(), N+1);
  EXPECT_DOUBLE_EQ(stationaryAcceptance.sum(), 1.0);
  for( unsigned int i=0; i<N+1; ++i ) { EXPECT_TRUE(stationaryAcceptance(i)>-1.0e-15); }

  // accept/reject M of the proposed states
  std::vector<std::shared_ptr<SamplingState> > newStates = kern->Step(0, state);
  EXPECT_EQ(newStates.size(), M);
}
