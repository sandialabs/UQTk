#include "MUQ/SamplingAlgorithms/MCMCFactory.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

std::shared_ptr<SingleChainMCMC>
MCMCFactory::CreateSingleChain(pt::ptree pt,
                               std::shared_ptr<AbstractSamplingProblem> problem) {
  std::shared_ptr<SingleChainMCMC> singleChain = std::make_shared<SingleChainMCMC>(pt, problem);

  return singleChain;
}
