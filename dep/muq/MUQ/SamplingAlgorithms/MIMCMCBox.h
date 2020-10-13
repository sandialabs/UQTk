#ifndef MIMCMCBOX_H
#define MIMCMCBOX_H

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/MIKernel.h"
#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <fstream>

namespace pt = boost::property_tree;
using namespace muq::Utilities;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Representation of a Multiindex MCMC telescoping sum component.
        @details This holds Markov chains whose differences are estimated
        in one component of the telescoping sum of an MIMCMC method.
        In the multilevel case, this contains two neighboring chains and
        allows computing their difference in expectation value; for multiindex,
        it extends to arbitrary dimensions.
     */
    class MIMCMCBox {
    public:

      MIMCMCBox(std::shared_ptr<MIComponentFactory> componentFactory, std::shared_ptr<MultiIndex> boxHighestIndex);

      void Sample();

      Eigen::VectorXd MeanQOI();

      void DrawChain(std::shared_ptr<SingleChainMCMC> chain, std::string chainid, std::ofstream& graphfile) const;

      void Draw(std::ofstream& graphfile, bool drawSamples = true) const;

      std::shared_ptr<SingleChainMCMC> FinestChain();

      std::shared_ptr<AbstractSamplingProblem> GetFinestProblem();

    private:

      // Creates a path of multiindices back to zero, preferring a route along the diagonal
      std::shared_ptr<MultiIndexSet> CreateRootPath(std::shared_ptr<MultiIndex> index);


      std::shared_ptr<MIComponentFactory> componentFactory;
      std::shared_ptr<MultiIndex> boxHighestIndex;
      std::shared_ptr<MultiIndex> boxLowestIndex;
      std::shared_ptr<MultiIndexSet> boxIndices;
      std::vector<std::shared_ptr<SingleChainMCMC>> boxChains;
      std::vector<std::shared_ptr<SingleChainMCMC>> tailChains;
      std::shared_ptr<AbstractSamplingProblem> finestProblem;
    };

  }
}

#endif
