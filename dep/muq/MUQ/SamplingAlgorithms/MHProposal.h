#ifndef MHPROPOSAL_H_
#define MHPROPOSAL_H_

#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class MHProposal
        @brief Implementation of the classic Random Walk Metropolis proposal
        @details <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "ProposalVariance"  | Double | - | The variance of an isotropic random walk proposal. |
    */
    class MHProposal : public MCMCProposal {
    public:

      MHProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob);

      MHProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob,
                 std::shared_ptr<muq::Modeling::GaussianBase> proposalIn);

      virtual ~MHProposal() = default;

    protected:

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::GaussianBase> proposal;

      virtual std::shared_ptr<SamplingState>
      Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double
      LogDensity(std::shared_ptr<SamplingState> currState,
                 std::shared_ptr<SamplingState> propState) override;


    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
