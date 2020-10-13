#include "MUQ/SamplingAlgorithms/MHProposal.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/AnyHelpers.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(MHProposal)

MHProposal::MHProposal(pt::ptree const& pt,
                       std::shared_ptr<AbstractSamplingProblem> prob) :
                       MCMCProposal(pt,prob) {

  unsigned int problemDim = prob->blockSizes(blockInd);

  // compute the (diagonal) covariance for the proposal
  const Eigen::VectorXd cov = pt.get("ProposalVariance", 1.0)*
                              Eigen::VectorXd::Ones(problemDim);

  // created a Gaussian with scaled identity covariance
  proposal = std::make_shared<Gaussian>(Eigen::VectorXd::Zero(problemDim), cov);
}

MHProposal::MHProposal(pt::ptree const& pt,
                       std::shared_ptr<AbstractSamplingProblem> prob,
                       std::shared_ptr<GaussianBase> proposalIn) :
                       MCMCProposal(pt,prob), proposal(proposalIn) {}

std::shared_ptr<SamplingState> MHProposal::Sample(std::shared_ptr<SamplingState> currentState) {
  assert(currentState->state.size()>blockInd);

  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;
  assert(props.size()>blockInd);
  Eigen::VectorXd const& xc = currentState->state.at(blockInd);

  Eigen::VectorXd prop = proposal->Sample();
  props.at(blockInd) = xc + prop;

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double MHProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                              std::shared_ptr<SamplingState> propState) {

  Eigen::VectorXd diff = propState->state.at(blockInd)-currState->state.at(blockInd);
  return proposal->LogDensity(diff);//, std::pair<boost::any, Gaussian::Mode>(conditioned->state.at(blockInd), Gaussian::Mode::Mean));
}
