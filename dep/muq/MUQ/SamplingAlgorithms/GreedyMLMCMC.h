#ifndef GreedyMLMCMC_H
#define GreedyMLMCMC_H

#include <boost/property_tree/ptree.hpp>

#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"
#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

namespace pt = boost::property_tree;

namespace muq {
  namespace SamplingAlgorithms {

    /** @brief Greedy Multilevel MCMC method.
        @details A Multilevel MCMC method choosing
        the number of samples adaptively at runtime,
        estimating the most profitable level from
        statistical information on samples.
     */
    class GreedyMLMCMC : public SamplingAlgorithm {
    public:
      GreedyMLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory);

      virtual std::shared_ptr<SampleCollection> GetSamples() const override;
      virtual std::shared_ptr<SampleCollection> GetQOIs() const override;

      Eigen::VectorXd MeanQOI();

      void Draw(bool drawSamples = true);

      std::shared_ptr<MIMCMCBox> GetBox(int index);

    protected:
      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) override;

    private:
      const double e;
      const double beta;
      const int levels;
      std::shared_ptr<MIComponentFactory> componentFactory;
      const int numInitialSamples;
      std::vector<std::shared_ptr<MIMCMCBox>> boxes;
      int verbosity;
    };

  }
}

#endif
