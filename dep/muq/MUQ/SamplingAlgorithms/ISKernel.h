#ifndef ISKERNEL_H_
#define ISKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

namespace muq {
  namespace SamplingAlgorithms {

    /// Importance sampling transition kernel
    /**
       Propose from a biasing distribution, compute the weight and return a state with those values
     */
    class ISKernel : public TransitionKernel {
    public:

      ISKernel(boost::property_tree::ptree const& pt, std::shared_ptr<SamplingProblem> problem);

      ~ISKernel();

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
    
  } // namespace SamplingAlgorithms
} // namespace muq


#endif
