#ifndef GMHKERNEL_H_
#define GMHKERNEL_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Queue.h>
#endif

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// A kernel for the generalized Metropolis-Hastings kernel
    /**
       Reference: "A general construction for parallelizing Metropolis-Hastings algorithms" (Calderhead, 2014)
     */
    class GMHKernel : public MHKernel {
    public:

      /**
	 @param[in] pt Options for this kenel and the standard Metropolis-Hastings kernel
	 @param[in] problem The problem we want to sample
       */
      GMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem);

      /**
	 @param[in] pt Options for this kenel and the standard Metropolis-Hastings kernel
	 @param[in] problem The problem we want to sample
	 @param[in] proposalIn The proposal for the MCMC chain
       */
      GMHKernel(boost::property_tree::ptree const& pt, std::shared_ptr<AbstractSamplingProblem> problem, std::shared_ptr<MCMCProposal> proposalIn);

      virtual ~GMHKernel();

      /**
	 Propose GMHKernel::N points and compute the cumulative distribution of the stationary distribution for the acceptance probability (GMHKernel::proposedStates and GMHKernel::stationaryAcceptance, respectively)
	 @param[in] t The current step in the MCMC chain
	 @param[in] state The current MCMC state
       */
      virtual void PreStep(unsigned int const t, std::shared_ptr<SamplingState> state) override;


      /**
	 Sample GMHKernel::M states from the proposed states GMHKernel::proposedStates.
	 @param[in] t The current step in the MCMC chain
	 @param[in] state The current MCMC state
	 \return The accepted states
       */
      virtual std::vector<std::shared_ptr<SamplingState> > Step(unsigned int const t, std::shared_ptr<SamplingState> state) override;

      /// Sample the stationary distribution
      /**
	 \return The new points
       */
      std::vector<std::shared_ptr<SamplingState> > SampleStationary() const;

      /// Get the cumulative stationary acceptance probability
      /**
	 Must be called after GMHKernel::PreStep.
	 \return The stationary acceptance probability
       */
      Eigen::VectorXd StationaryAcceptance() const;

#if MUQ_HAS_PARCER
      std::shared_ptr<parcer::Communicator> GetCommunicator() const;
#endif

    private:

      /// Propose \f$N\f$ points in serial and evaluate the log target
      /**
	 @param[in] t The current step in the MCMC chain
	 @param[in] state The current point
       */
      void SerialProposal(unsigned int const t, std::shared_ptr<SamplingState> state);

#if MUQ_HAS_PARCER
      /// Propose \f$N\f$ points in parallel and evaluate the log target
      /**
	 @param[in] t The current step in the MCMC chain
	 @param[in] state The current point
       */
      void ParallelProposal(unsigned int const t, std::shared_ptr<SamplingState> state);
#endif

      /// Compute the stationary transition density
      /**
	 @param[in] R The log-target
       */
      void AcceptanceDensity(Eigen::VectorXd& R);
      
      Eigen::MatrixXd AcceptanceMatrix(Eigen::VectorXd const& R) const;

      /// Compute the cumulative acceptance density
      /**
	 @param[in] R The log-target plus the log-proposals
       */
      void ComputeStationaryAcceptance(Eigen::VectorXd const& R);

      /// Number of proposals
      const unsigned int N;

      /// Number of proposals plus one
      /**
	 Defined so we don't aways have to compute \f$N+1\f$.
       */
      const unsigned int Np1;

      /// Number of accepted points (number of points added to the chain)
      const unsigned int M;

      /// The cumulative stationary accepatnce probability 
      Eigen::VectorXd stationaryAcceptance;

      /// Proposed states
      std::vector<std::shared_ptr<SamplingState> > proposedStates;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
