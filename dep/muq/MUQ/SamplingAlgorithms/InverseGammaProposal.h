#ifndef INVERSEGAMMAPROPOSAL_H_
#define INVERSEGAMMAPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
      @ingroup MCMCProposals
      @class InverseGammaProposal
      @brief Defines a proposal using the analytic conditional Inverse Gamma distribution for the variance of a Gaussian distribution
      @details Consider a Metropolis-Within-Gibbs sampler for a problem where
      the inverse Gamma distribution is used to model the variance of a Gaussian
      distribution.  In that setting, the distribution of the variance given the
      state of the Gaussian random variable is known analytically and can be sampled
      from directly.  This proposal leverages that fact. It assumes a prior distribution
      over the variance is given by \f$\sigma^2 \sim IG(\alpha,\beta)\f$ and that
      the Gaussian random variable has zero mean and covariance \f$\simga^2 I\f$.  Then,
      given an observation of the Gaussian random variable \f$x=[x_1,x_2, \ldots, x_N]\f$,
      we know that
      \f[
      \sigma^2 | x \sim IG(\alpha + \frac{N}{2}, \beta + \frac{1}{2}\sum_{i=1}^Nx_i^2.
      \f]
      This is the proposal density defined by this class.
    */
    class InverseGammaProposal : public MCMCProposal {
    public:

      InverseGammaProposal(boost::property_tree::ptree       const& pt,
                           std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~InverseGammaProposal() = default;

    protected:

      /// The prior value of alpha
      const double alpha;

      /// The prior value of beta
      const double beta;

      /// The index of the Gaussian block
      const unsigned int gaussBlock;

      /// The mean of the Gaussian distribution
      const Eigen::VectorXd gaussMean;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override;


      static Eigen::VectorXd ExtractMean(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& gaussNode);

      static std::shared_ptr<muq::Modeling::InverseGamma> ExtractInverseGamma(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static double ExtractAlpha(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static double ExtractBeta(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& igNode);
      static unsigned int ExtractGaussInd(std::shared_ptr<AbstractSamplingProblem> prob, std::string const& gaussNode);
    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
