#include "MUQ/Utilities/AnyHelpers.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;


REGISTER_MCMC_PROPOSAL(AMProposal)
AMProposal::AMProposal(pt::ptree const& pt , std::shared_ptr<AbstractSamplingProblem> prob) : MHProposal(pt, prob),
                                              adaptSteps(pt.get<unsigned int>("AdaptSteps")),
                                              adaptStart(pt.get<unsigned int>("AdaptStart")),
                                              adaptScale(pt.get<double>("AdaptScale")) {}


AMProposal::~AMProposal() {}

void AMProposal::Adapt(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& states) {
  // always update the sample mean and covariance
  Update(t, states);

  if( t%adaptSteps==0 && t>adaptStart ) {

    Eigen::MatrixXd adjustedCov = adaptScale * cov + 1e-10 * Eigen::MatrixXd::Identity(cov.rows(), cov.cols());

    // update the proposal covariance
    std::dynamic_pointer_cast<Gaussian>(proposal)->SetCovariance(adjustedCov);
  }
}

void AMProposal::UpdateOne(unsigned int const numSamps, std::shared_ptr<SamplingState> state)
{
  // first sample---we have no mean, just set it to the first sample
  if( mean.size()==0 ){
    mean = state->state.at(blockInd);
    return;
  }

  // update the mean
  Eigen::VectorXd oldMean = mean;
  Eigen::VectorXd const& newState  = state->state.at(blockInd);

  mean = (oldMean*numSamps + newState)/(numSamps+1.0);

  // If we haven't compute the covariance before...
  if( cov.rows()==0 ){

    cov = Eigen::MatrixXd::Zero(oldMean.size(),oldMean.size());

    //compute covariance from scratch, from the definition
    cov.selfadjointView<Eigen::Lower>().rankUpdate(oldMean - mean, 1.0);
    cov.selfadjointView<Eigen::Lower>().rankUpdate(newState - mean, 1.0);

  }else{
    cov *= (numSamps - 1.0) / numSamps;
    //note that the asymmetric form fixes the fact that the old mean was wrong
    cov += (1.0 / static_cast<double>(numSamps)) * (newState - oldMean) * (newState - mean).transpose();
  }
}

void AMProposal::Update(unsigned int const numSamps, std::vector<std::shared_ptr<SamplingState>> const& states) {

  for(int i=0; i<states.size(); ++i)
    UpdateOne(numSamps+i,states.at(i));

}
