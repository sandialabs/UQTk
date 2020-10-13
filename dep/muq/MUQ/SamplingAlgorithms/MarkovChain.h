#ifndef MarkovChain_H
#define MarkovChain_H

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms{

    /**
    @ingroup MCMC
    @class MarkovChain
    @brief A class for storing and working with the results of Markov chain Monte Carlo algorithms.
    @details The MarkovChain class is a child of SampleCollection where the sample
    weights correspond to the number of consecutive steps taking the same value,
    and the weights are unnormalized (i.e., do not sum to one).  This is a useful
    class for storing the chain produced by an MCMC algorithm without storing the
    duplicate points that result from rejected proposals.
    */
    class MarkovChain : public SampleCollection
    {
    public:

      MarkovChain() = default;

      virtual ~MarkovChain() = default;

      /** Computes the effective sample size using the method described
      Ulli Wolff's "Monte Carlo errors with less error"
      */
      virtual Eigen::VectorXd ESS(int blockDim=-1) const override;

      static double SingleComponentESS(Eigen::Ref<const Eigen::VectorXd> const& trace);

    private:

      std::vector<std::unordered_map<std::string, boost::any> > meta;

    }; // class MarkovChain
  }
}

#endif
