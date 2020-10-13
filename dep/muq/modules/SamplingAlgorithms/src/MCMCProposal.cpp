#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

MCMCProposal::MCMCProposal(boost::property_tree::ptree       const& pt,
                           std::shared_ptr<AbstractSamplingProblem> prob) : blockInd(pt.get("BlockIndex",0)) {}

std::shared_ptr<MCMCProposal> MCMCProposal::Construct(pt::ptree                         const& pt,
                                                      std::shared_ptr<AbstractSamplingProblem> prob) {

  // get the name of the proposal
  std::string proposalName = pt.get<std::string>("Method");

  auto proposalMap = GetMCMCProposalMap();
  auto iter = proposalMap->find(proposalName);

  if(iter == proposalMap->end()){
    std::cerr << "ERROR: Could not find MCMC proposal \"" << proposalName << "\".  Available options are:\n";

    for(auto it=proposalMap->begin(); it!=proposalMap->end(); ++it)
      std::cerr << "  " << it->first << std::endl;
    std::cerr << std::endl;

    assert(iter!=proposalMap->end());
  }

  return iter->second(pt,prob);

}

std::shared_ptr<MCMCProposal::MCMCProposalMap> MCMCProposal::GetMCMCProposalMap() {
  // define a static map from type to constructor
  static std::shared_ptr<MCMCProposalMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<MCMCProposalMap>();
  }

  return map;
}

#if MUQ_HAS_PARCER
void MCMCProposal::SetCommunicator(std::shared_ptr<parcer::Communicator> newcomm) {
  comm = newcomm;
}
#endif
