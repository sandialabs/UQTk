#include "MUQ/SamplingAlgorithms/ParallelAMProposal.h"

#if MUQ_HAS_PARCER

namespace pt = boost::property_tree;
using namespace muq::SamplingAlgorithms;

REGISTER_MCMC_PROPOSAL(ParallelAMProposal)
ParallelAMProposal::ParallelAMProposal(pt::ptree const& pt , std::shared_ptr<AbstractSamplingProblem> prob) : AMProposal(pt, prob) {}

void ParallelAMProposal::Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState> > const& states) {
  assert(comm);

  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    // get the samples from the other processors (or send your own to them)
    std::vector<std::shared_ptr<SamplingState> > otherStates;
    if( i==comm->GetRank() ) { otherStates = states; }
    comm->Bcast(otherStates, i);

    // adapt
    totSamps += otherStates.size();
    AMProposal::Adapt(totSamps, otherStates);
  }
}

#endif // end MUQ_HAS_PARCER
