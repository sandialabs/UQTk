#ifndef MCMCPROPOSAL_H_
#define MCMCPROPOSAL_H_

#include <map>

#include <functional>
#include <boost/property_tree/ptree.hpp>

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Communicator.h>
#endif

#include "MUQ/Utilities/RegisterClassName.h"

#include "MUQ/Modeling/Distributions/Distribution.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"
#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    class MCMCProposal : public std::enable_shared_from_this<MCMCProposal> {
    public:

      MCMCProposal(boost::property_tree::ptree       const& pt,
                   std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~MCMCProposal() = default;

      /// Static constructor for the transition kernel
      /**
	 @param[in] pt The options for the MCMC kernel
	 \return The MCMC proposal
       */
      static std::shared_ptr<MCMCProposal> Construct(boost::property_tree::ptree       const& pt,
                                                     std::shared_ptr<AbstractSamplingProblem> prob);

      typedef std::function<std::shared_ptr<MCMCProposal>(boost::property_tree::ptree,std::shared_ptr<AbstractSamplingProblem>)> MCMCProposalConstructor;
      typedef std::map<std::string, MCMCProposalConstructor> MCMCProposalMap;
      static std::shared_ptr<MCMCProposalMap> GetMCMCProposalMap();

      /// Adapt the proposal after each step
      /**
	 By default this function does nothing but children can override it to adapt the proposal
	 @param[in] t The current step
	 @param[in] state The current state
       */
      virtual void Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state) {};

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) = 0;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) = 0;

#if MUQ_HAS_PARCER
      void SetCommunicator(std::shared_ptr<parcer::Communicator> newcomm);
#endif

      const int blockInd = 0;

    protected:
      
#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> comm;
#endif
    };
  } // namespace SamplingAlgoirthms
} // namespace muq

#define REGISTER_MCMC_PROPOSAL(NAME) static auto reg ##NAME		\
  = muq::SamplingAlgorithms::MCMCProposal::GetMCMCProposalMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));


#endif
