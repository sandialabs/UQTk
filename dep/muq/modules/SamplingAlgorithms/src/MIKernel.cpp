#include "MUQ/SamplingAlgorithms/MIKernel.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

MIKernel::MIKernel(boost::property_tree::ptree const& pt,
                   std::shared_ptr<AbstractSamplingProblem> problem,
                   std::shared_ptr<AbstractSamplingProblem> coarse_problem,
                   std::shared_ptr<MCMCProposal> proposal,
                   std::shared_ptr<MCMCProposal> coarse_proposal,
                   std::shared_ptr<MIInterpolation> proposalInterpolation,
                   std::shared_ptr<SingleChainMCMC> coarse_chain)
  : TransitionKernel(pt, problem),
    coarse_problem(coarse_problem),
    proposal(proposal),
    coarse_proposal(coarse_proposal),
    proposalInterpolation(proposalInterpolation),
    coarse_chain(coarse_chain) {}


MIKernel::~MIKernel() {}

void MIKernel::PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state){
  proposal->Adapt(t,state);
}

std::vector<std::shared_ptr<SamplingState>> MIKernel::Step(unsigned int const t, std::shared_ptr<SamplingState> prevState){
  assert(proposal);

  // If no coarse sample is specified, assume it's the first one
  std::shared_ptr<SamplingState> coarsePrevState;
  if(prevState->HasMeta("coarseSample")) {
    coarsePrevState = AnyCast(prevState->meta["coarseSample"]);
  } else {
    coarsePrevState = coarse_chain->GetSamples()->at(0);
  }

  // New fine proposal
  std::shared_ptr<SamplingState> prop = proposal->Sample(prevState);
  std::shared_ptr<SamplingState> coarseProp = coarse_proposal->Sample(coarsePrevState);
  std::shared_ptr<SamplingState> fineProp = proposalInterpolation->Interpolate (coarseProp, prop);

  // compute acceptance probability
  double propTarget;
  double currentTarget;
  double propTargetCoarse;
  double currentTargetCoarse;

  if(prevState->HasMeta("LogTarget")){
    currentTarget = AnyCast(prevState->meta["LogTarget"]);
  }else{
    currentTarget = problem->LogDensity(t, prevState, AbstractSamplingProblem::SampleType::Accepted);
    prevState->meta["LogTarget"] = currentTarget;
    if (problem->numBlocksQOI > 0) {
      prevState->meta["QOI"] = problem->QOI();
    }
  }

  if(coarsePrevState->HasMeta("LogTarget")){
    currentTargetCoarse = AnyCast(coarsePrevState->meta["LogTarget"]);
  }else{
    currentTargetCoarse = coarse_problem->LogDensity(t, coarsePrevState, AbstractSamplingProblem::SampleType::Accepted);
    coarsePrevState->meta["LogTarget"] = currentTargetCoarse;
  }

  propTarget = problem->LogDensity(t, fineProp, AbstractSamplingProblem::SampleType::Proposed);
  if(coarseProp->HasMeta("LogTarget")){
    propTargetCoarse = AnyCast(coarseProp->meta["LogTarget"]);
  }else{
    propTargetCoarse = coarse_problem->LogDensity(t, coarseProp, AbstractSamplingProblem::SampleType::Proposed);
    coarseProp->meta["LogTarget"] = propTargetCoarse;
  }
  fineProp->meta["LogTarget"] = propTarget;
  fineProp->meta["coarseSample"] = coarseProp; // Hook up fine sample to coarse one

  // Aceptance probability
  const double forwardPropDens = proposal->LogDensity(prevState, fineProp);
  const double backPropDens = proposal->LogDensity(fineProp, prevState);
  const double forwardPropDensCoarse = coarse_proposal->LogDensity(coarsePrevState, coarseProp);
  const double backPropDensCoarse = coarse_proposal->LogDensity(coarseProp, coarsePrevState);
  const double alpha = std::exp(propTarget + backPropDens - forwardPropDens - currentTarget
                              -(propTargetCoarse + backPropDensCoarse - forwardPropDensCoarse - currentTargetCoarse));

  // accept/reject
  numCalls++;
  if(RandomGenerator::GetUniform()<alpha) {
    numAccepts++;
    if (problem->numBlocksQOI > 0) {
      fineProp->meta["QOI"] = problem->QOI();
    }
    return std::vector<std::shared_ptr<SamplingState>>(1,fineProp);
  } else {
    // Return copy of previous state in order to attach the new coarse proposal to it
    auto prevStateCopy = std::make_shared<SamplingState>(prevState->state);
    prevStateCopy->meta["LogTarget"] = prevState->meta["LogTarget"];
    if (prevState->HasMeta("QOI")) {
      prevStateCopy->meta["QOI"] = prevState->meta["QOI"];
    }
    prevStateCopy->meta["coarseSample"] = coarseProp;
    return std::vector<std::shared_ptr<SamplingState>>(1,prevStateCopy);
  }
}

void MIKernel::PrintStatus(const std::string prefix) const
{
  std::stringstream msg;
  msg << std::setprecision(2);
  msg << prefix << "Acceptance Rate = "  << 100.0*double(numAccepts)/double(numCalls) << "%";

  std::cout << msg.str() << std::endl;
}
