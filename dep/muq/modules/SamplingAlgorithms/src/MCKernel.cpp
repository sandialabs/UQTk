#include "MUQ/SamplingAlgorithms/MCKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(MCKernel)

MCKernel::MCKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem) {}

std::vector<std::shared_ptr<SamplingState>>  MCKernel::Step(std::shared_ptr<SamplingState> prevState)
{
  const boost::any newState = problem->GetDistribution()->Sample(prevState->state);
  return std::vector<std::shared_ptr<SamplingState>>(1,std::make_shared<SamplingState>(newState, 1.0));
}
