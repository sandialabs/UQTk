#ifndef MCKERNEL_H_
#define MCKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace muq {
  namespace SamplingAlgorithms {
    /// Monte Carlo transition kernel
    /**
       Samples from the target distirbution directly and returns that state.
     */
    class MCKernel : public TransitionKernel {
    public:

      MCKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      virtual ~MCKernel() = default;

      virtual std::vector<std::shared_ptr<SamplingState>>  Step(std::shared_ptr<SamplingState> prevState) override;

    private:

      /// The number of Monte Carlo samples
      //const unsigned int N;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
