#ifndef MIINTERPOLATION_H_
#define MIINTERPOLATION_H_

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Interpolation interface combining coarse and fine samples.
        @details This interface defines how samples drawn from coarse chains
        and fine proposals are combined to form a fine chain proposal,
        as needed for Multiindex MCMC.
        In many cases, this simply means concatenating the associated state vectors.
     */
    class MIInterpolation {
    public:

      virtual ~MIInterpolation() = default;

      virtual std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) = 0;


    };


  }
}

#endif
