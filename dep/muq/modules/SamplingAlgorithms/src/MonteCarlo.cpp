#include "MUQ/SamplingAlgorithms/MonteCarlo.h"

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

MonteCarlo::MonteCarlo() :
  SamplingAlgorithm(false) // the sampling algorithm produces uncorrelated samples
{}

MonteCarlo::~MonteCarlo() {}

std::shared_ptr<TransitionKernel> MonteCarlo::ConstructKernel(pt::ptree& pt, std::shared_ptr<SamplingProblem> problem) const {
  // we have to use the Monte Carlo transition kernel
  pt.put<std::string>("SamplingAlgorithm.TransitionKernel", "MCKernel");

  return TransitionKernel::Construct(pt, problem);
}
