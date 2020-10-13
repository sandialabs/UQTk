#ifndef SINGLECHAINMCMC_H
#define SINGLECHAINMCMC_H

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"
#include "MUQ/SamplingAlgorithms/TransitionKernel.h"
#include "MUQ/SamplingAlgorithms/ThinScheduler.h"

#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace muq{
  namespace SamplingAlgorithms{

    /** @ingroup MCMC
        @class SingleChainMCMC
        @brief Defines an MCMC sampler with a single chain
        @details
        <B>Configuration Parameters:</B>
        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "NumSamples"  | Int | - | The total number of steps (including burnin) to take, i.e., the length of the Markov chain. |
        "BurnIn"      | Int | 0 | The number of steps at the beginning of the chain to ignore. |
        "PrintLevel"  | Int | 3 | The amount of information to print to std::cout. Valid values are in [0,1,2,3] with  0 = Nothing, 3 = The most |
        "KernelList"  | String | - | A comma separated list of other parameter blocks that define the transition kernels for each Metropolis-in-Gibbs block. |
    */
    class SingleChainMCMC : public SamplingAlgorithm
    {

    public:

#if MUQ_HAS_PARCER
      SingleChainMCMC(boost::property_tree::ptree pt,
                      std::shared_ptr<parcer::Communicator> const& comm,
                      std::vector<std::shared_ptr<TransitionKernel> > const& kernelsIn);

      SingleChainMCMC(boost::property_tree::ptree pt,
                      std::shared_ptr<AbstractSamplingProblem> const& problem,
                      std::shared_ptr<parcer::Communicator> const& comm);
#endif

      SingleChainMCMC(boost::property_tree::ptree pt,
                      std::shared_ptr<AbstractSamplingProblem> const& problem);

      SingleChainMCMC(boost::property_tree::ptree pt,
                      std::vector<std::shared_ptr<TransitionKernel> > const& kernels);

      virtual ~SingleChainMCMC() = default;

      /// Set the state of the MCMC chain
      /**
        If no steps have been taken, this function sets the starting point.
        Set the current state AND adds it to the sample collection.
      */
      virtual void SetState(std::vector<Eigen::VectorXd> const& x0) override;

      template<typename... Args>
        inline void SetState(Args const&... args) {
          std::vector<Eigen::VectorXd> vec;
          SetStateRecurse(vec, args...);
      }

      virtual std::vector<std::shared_ptr<TransitionKernel>>& Kernels(){return kernels;};

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

      virtual void Sample();

      virtual double TotalTime() { return totalTime; }

    protected:

      std::shared_ptr<SamplingState> SaveSamples(std::vector<std::shared_ptr<SamplingState> > const& newStates, std::shared_ptr<SamplingState>& lastSavedState, unsigned int& sampNum) const;

      bool ShouldSave(unsigned int const sampNum) const;

      void PrintStatus(unsigned int currInd) const{PrintStatus("",currInd);};
      void PrintStatus(std::string prefix, unsigned int currInd) const;

      std::shared_ptr<SaveSchedulerBase> scheduler;
      std::shared_ptr<SaveSchedulerBase> schedulerQOI;

      unsigned int numSamps;
      unsigned int burnIn;
      unsigned int printLevel;

      // A vector of transition kernels: One for each block
      std::vector<std::shared_ptr<TransitionKernel>> kernels;

    private:

      template<typename... Args>
        inline void SetStateRecurse(std::vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& it, Args const&... args) {
          vec.push_back(it);
          SetStateRecurse(vec, args...);
        }

      inline void SetStateRecurse(std::vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
        vec.push_back(last);
        SetState(vec);
      }

      unsigned int sampNum = 1;
      std::shared_ptr<SamplingState> prevState = nullptr;
      std::shared_ptr<SamplingState> lastSavedState = nullptr;
      double totalTime = 0.0;

      void Setup(boost::property_tree::ptree pt,
                 std::vector<std::shared_ptr<TransitionKernel>> const& kernelsIn);

      void Setup(boost::property_tree::ptree pt, std::shared_ptr<AbstractSamplingProblem> const& problem);
    }; // class SingleChainMCMC
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // #ifndef SINGLECHAINMCMC_H
