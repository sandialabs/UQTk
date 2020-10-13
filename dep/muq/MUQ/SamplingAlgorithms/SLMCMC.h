#ifndef SLMCMC_H
#define SLMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Single-level MCMC for multiindex sampling problems.
        @details A wrapper generating a single-chain MCMC
        based on the finest problem of a multiindex sampling problem.
        This is mostly needed for computing reference solutions to
        multilevel/multiindex MCMC methods.
    */

    class SLMCMC : public SamplingAlgorithm {

    public:
      SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      virtual std::shared_ptr<SampleCollection> GetSamples() const override;
      virtual std::shared_ptr<SampleCollection> GetQOIs() const override;

      Eigen::VectorXd MeanQOI();

      Eigen::VectorXd MeanParameter();

    protected:
      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

    private:
      std::shared_ptr<MIComponentFactory> componentFactory;

      std::shared_ptr<SingleChainMCMC> coarse_chain;
    };

  }
}

#endif
