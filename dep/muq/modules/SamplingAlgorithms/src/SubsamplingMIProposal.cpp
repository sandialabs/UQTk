#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"
namespace muq {
  namespace SamplingAlgorithms {

    SubsamplingMIProposal::SubsamplingMIProposal (pt::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<SingleChainMCMC> coarseChain)
     : MCMCProposal(pt,prob), coarseChain(coarseChain),
       subsampling(pt.get("subsampling",1))
    {}

    std::shared_ptr<SamplingState> SubsamplingMIProposal::Sample(std::shared_ptr<SamplingState> currentState) {

      sampleID += subsampling;
      while (coarseChain->GetSamples()->size() <= sampleID) {
        coarseChain->Sample();
      }

      return coarseChain->GetSamples()->at(sampleID);
    }

    double SubsamplingMIProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                                             std::shared_ptr<SamplingState> propState) {
      return 0;
    }


  }
}
