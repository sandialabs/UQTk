#ifndef MICOMPONENTFACTORY_H
#define MICOMPONENTFACTORY_H

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MIInterpolation.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

using namespace muq::Utilities;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Interface defining models on a multiindex structure.
        @details Defines all components of a multiindex model. MLMCMC and MIMCMC
        work on this structure.
        In particular, a sampling problem is defined for each multiindex along with
        associated proposal densities and other required components.
     */
    class MIComponentFactory {
    public:
      virtual ~MIComponentFactory() = default;
      virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) = 0;
      virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const&index,
                                                            std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                            std::shared_ptr<SingleChainMCMC> const& coarseChain) = 0;
      virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) = 0;
      virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) = 0;
      virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) = 0;
      virtual std::shared_ptr<MultiIndex> FinestIndex() = 0;

    };

  }
}

#endif
