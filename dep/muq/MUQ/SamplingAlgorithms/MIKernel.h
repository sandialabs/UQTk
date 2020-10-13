#ifndef MIKERNEL_H_
#define MIKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

#include "MUQ/SamplingAlgorithms/MIInterpolation.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief MCMC kernel for Multiindex methods.
        @details This kernel combines a coarse proposal from a coarse chain
        with a fine one, as needed for MIMCMC.
     */
    class MIKernel : public TransitionKernel {
    public:

      MIKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem,
               std::shared_ptr<AbstractSamplingProblem> coarse_problem,
               std::shared_ptr<MCMCProposal> proposal,
               std::shared_ptr<MCMCProposal> coarse_proposal,
               std::shared_ptr<MIInterpolation> proposalInterpolation,
               std::shared_ptr<SingleChainMCMC> coarse_chain);

      ~MIKernel();

      virtual std::shared_ptr<MCMCProposal> Proposal(){return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual void PrintStatus(std::string prefix) const override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState) override;

      virtual double AcceptanceRate() const{return double(numAccepts)/double(numCalls);};

    protected:
      std::shared_ptr<AbstractSamplingProblem> coarse_problem;
      std::shared_ptr<MCMCProposal> proposal;
      std::shared_ptr<MCMCProposal> coarse_proposal;
      std::shared_ptr<MIInterpolation> proposalInterpolation;
      std::shared_ptr<SingleChainMCMC> coarse_chain;

      unsigned int numCalls = 0;
      unsigned int numAccepts = 0;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
