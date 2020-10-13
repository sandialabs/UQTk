#include "MUQ/SamplingAlgorithms/MCMC.h"

using namespace muq::SamplingAlgorithms;

MCMC::MCMC() :
  SamplingAlgorithm(true) // the sampling algorithm produces correlated samples
{}

MCMC::~MCMC() {}

std::shared_ptr<TransitionKernel> MCMC::ConstructKernel(boost::property_tree::ptree& pt, std::shared_ptr<SamplingProblem> problem) const {
  return TransitionKernel::Construct(pt, problem);
}
