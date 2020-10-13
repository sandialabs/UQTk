#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace muq {
  namespace SamplingAlgorithms {

    class TransitionKernel;

    class MonteCarlo : public SamplingAlgorithm {
    public:

      MonteCarlo();

      ~MonteCarlo();

    private:

      /// Create the transition kernel
      /**
	 @param[in] pt Parameters for the kernel
	 @param[in] problem The sampling problem that samples the distribution we are trying to characterize
	 \return The transition kernel
       */
      virtual std::shared_ptr<TransitionKernel> ConstructKernel(boost::property_tree::ptree&             pt,
                                                                std::shared_ptr<AbstractSamplingProblem> problem) const override;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
