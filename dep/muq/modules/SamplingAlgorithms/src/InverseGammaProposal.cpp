#include "MUQ/SamplingAlgorithms/InverseGammaProposal.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/WorkGraph.h"

#include "MUQ/Utilities/Exceptions.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

REGISTER_MCMC_PROPOSAL(InverseGammaProposal)

InverseGammaProposal::InverseGammaProposal(pt::ptree                         const& pt,
                                           std::shared_ptr<AbstractSamplingProblem> prob) : MCMCProposal(pt,prob),
                                                                                            alpha(ExtractAlpha(prob, pt.get<std::string>("InverseGammaNode"))),
                                                                                            beta(ExtractBeta(prob, pt.get<std::string>("InverseGammaNode"))),
                                                                                            gaussBlock(ExtractGaussInd(prob, pt.get<std::string>("GaussianNode"))),
                                                                                            gaussMean(ExtractMean(prob, pt.get<std::string>("GaussianNode"))) {
  unsigned int problemDim = prob->blockSizes(blockInd);
  if(problemDim != 1)
    throw muq::WrongSizeError("The InverseGammaProposal is only defined for a block of size 1.");
}

std::shared_ptr<SamplingState> InverseGammaProposal::Sample(std::shared_ptr<SamplingState> currentState) {


  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;
  Eigen::VectorXd const& gaussState = currentState->state.at(gaussBlock);

  double newAlpha = alpha + 0.5*gaussState.size();
  double newBeta = beta + 0.5*(gaussState - gaussMean).sum();
  props.at(blockInd) = InverseGamma(newAlpha, newBeta).Sample();

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double InverseGammaProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                                        std::shared_ptr<SamplingState> propState) {

  Eigen::VectorXd const& gaussState = currState->state.at(gaussBlock);
  Eigen::VectorXd const& sigmaState = propState->state.at(blockInd);

  double newAlpha = alpha + 0.5*gaussState.size();
  double newBeta = beta + 0.5*(gaussState - gaussMean).sum();

  return InverseGamma(newAlpha, newBeta).LogDensity(sigmaState);
}

Eigen::VectorXd InverseGammaProposal::ExtractMean(std::shared_ptr<AbstractSamplingProblem> prob,
                                                  std::string                        const& gaussNode)
{
  // Cast the abstract base class into a sampling problem
  std::shared_ptr<SamplingProblem> prob2 = std::dynamic_pointer_cast<SamplingProblem>(prob);
  assert(prob2);

  // From the sampling problem, extract the ModPiece and try to cast it to a ModGraphPiece
  std::shared_ptr<ModPiece> targetDens = prob2->GetDistribution();
  std::shared_ptr<ModGraphPiece> targetDens2 = std::dynamic_pointer_cast<ModGraphPiece>(targetDens);
  assert(targetDens2);

  // Get the graph
  auto graph = targetDens2->GetGraph();

  // Get the Gaussian node
  auto gaussPiece = graph->GetPiece(gaussNode);

  // try to cast the gaussPiece to a density
  auto dens = std::dynamic_pointer_cast<Density>(gaussPiece);
  assert(dens);

  auto gaussDist = std::dynamic_pointer_cast<Gaussian>(dens->GetDistribution());
  assert(gaussDist);

  return gaussDist->GetMean();
}

std::shared_ptr<InverseGamma> InverseGammaProposal::ExtractInverseGamma(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode)
{
  // Cast the abstract base class into a sampling problem
  std::shared_ptr<SamplingProblem> prob2 = std::dynamic_pointer_cast<SamplingProblem>(prob);
  assert(prob2);

  // From the sampling problem, extract the ModPiece and try to cast it to a ModGraphPiece
  std::shared_ptr<ModPiece> targetDens = prob2->GetDistribution();
  std::shared_ptr<ModGraphPiece> targetDens2 = std::dynamic_pointer_cast<ModGraphPiece>(targetDens);
  assert(targetDens2);

  // Get the graph
  auto graph = targetDens2->GetGraph();

  // Get the Gaussian node
  auto piece = graph->GetPiece(igNode);

  // try to cast the gaussPiece to a density
  auto dens = std::dynamic_pointer_cast<Density>(piece);
  assert(dens);

  auto dist = std::dynamic_pointer_cast<InverseGamma>(dens->GetDistribution());
  assert(dist);

  return dist;
}

double InverseGammaProposal::ExtractAlpha(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode)
{
  return ExtractInverseGamma(prob,igNode)->alpha;
}
double InverseGammaProposal::ExtractBeta(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode)
{
  return ExtractInverseGamma(prob,igNode)->beta;
}
unsigned int InverseGammaProposal::ExtractGaussInd(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& gaussNode)
{
  // Cast the abstract base class into a sampling problem
  std::shared_ptr<SamplingProblem> prob2 = std::dynamic_pointer_cast<SamplingProblem>(prob);
  assert(prob2);

  // From the sampling problem, extract the ModPiece and try to cast it to a ModGraphPiece
  std::shared_ptr<ModPiece> targetDens = prob2->GetDistribution();
  std::shared_ptr<ModGraphPiece> targetDens2 = std::dynamic_pointer_cast<ModGraphPiece>(targetDens);
  assert(targetDens2);

  // Get the graph
  auto graph = targetDens2->GetGraph();

  // Get the piece corresponding to the Gaussian name
  auto gaussPiece = graph->GetPiece(gaussNode + "_0");
  assert(gaussPiece);

  // Figure out which one of the input pieces corresponds to the Gaussian
  std::vector<std::shared_ptr<ConstantVector> > constantPieces = targetDens2->GetConstantPieces();

  auto iter = std::find_if(constantPieces.begin(), constantPieces.end(), [gaussPiece](std::shared_ptr<ModPiece> piece)->bool{return gaussPiece==piece;} );
  assert(iter!=constantPieces.end());

  return std::distance(iter, constantPieces.begin());
}
