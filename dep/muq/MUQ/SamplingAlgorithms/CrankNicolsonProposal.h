#ifndef CRANKNICOLSONPROPOSAL_H_
#define CRANKNICOLSONPROPOSAL_H_

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /**
    @ingroup MCMCProposals
    @class CrankNicolsonProposal
    @brief An implement of the dimension-independent pCN proposal.

    @details
    This class implements the preconditioned Crank Nicolson proposal (pCN) from
    Cotter et al., 2013.  The proposal takes the form
    \f[
    u^\prime = \sqrt{ 1 - \beta^2} u_c + \beta z,
    \f]
    where \f$u_c\f$ is the current state of the chain, \f$z\sim N(0,C)\f$ is a normal
    random variable with a strategically chosen covariance \f$C\f$ (often the prior covariance), and \f$u^\prime\f$
    is the propsed point.  The parameter \f$\beta\f$ is a tuning parameter.

    <B>Configuration Parameters:</B>
    Parameter Key | Type | Default Value | Description |
    ------------- | ------------- | ------------- | ------------- |
    "Beta"  | Double | 0.5 | The proposal scaling \f$\beta\f$ defined above. |
    "PriorNode" | String | - | (Optional.)  If specified, this class assumes the target density was constructed from a WorkGraph and will look set the value of the covariance \f$C\f$ to the covariance of the Gaussian density at the specified node. |
    */
    class CrankNicolsonProposal : public MCMCProposal {
    public:

      CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                            std::shared_ptr<AbstractSamplingProblem> prob,
                            std::shared_ptr<muq::Modeling::GaussianBase> prior);

      CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                            std::shared_ptr<AbstractSamplingProblem> prob);

      virtual ~CrankNicolsonProposal() = default;

    protected:

      double beta;

      // When the prior distribution has inputs, we need to save information to evaluate them
      std::shared_ptr<muq::Modeling::ModPiece> priorMeanModel;
      std::vector<int> priorMeanInds;

      std::shared_ptr<muq::Modeling::ModPiece> priorCovModel;
      std::vector<int> priorCovInds;
      bool priorUsesCov;

      std::vector<Eigen::VectorXd>GetPriorInputs(std::vector<Eigen::VectorXd> const& currState);

      //const Eigen::VectorXd priorMu;

      /// The proposal distribution
      //std::shared_ptr<muq::Modeling::Gaussian> propPart;

      // Sometimes we have to keep track of the prior distribution so we can update the proposal mean and covariance
      std::shared_ptr<muq::Modeling::GaussianBase> priorDist;

      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> currentState) override;

      virtual double LogDensity(std::shared_ptr<SamplingState> currState,
                                std::shared_ptr<SamplingState> propState) override;

      void ExtractPrior(std::shared_ptr<AbstractSamplingProblem> prob, std::string nodeName);
    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
