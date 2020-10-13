#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

GaussianBase::GaussianBase(unsigned int dim) : Distribution(dim, Eigen::VectorXi()),
                                               mean(Eigen::VectorXd::Zero(dim))
{}

GaussianBase::GaussianBase(unsigned int dim,
                           Eigen::VectorXi const& hyperSizesIn) : Distribution(dim, hyperSizesIn),
                                                                  mean(Eigen::VectorXd::Zero(dim))
{}



GaussianBase::GaussianBase(Eigen::VectorXd const& muIn) : Distribution(muIn.size(), Eigen::VectorXi()),
                                                         mean(muIn)
{}


GaussianBase::GaussianBase(Eigen::VectorXd const& muIn,
                   Eigen::VectorXi const& hyperSizesIn) : Distribution(muIn.size(), hyperSizesIn),
                                                         mean(muIn)
{}


double GaussianBase::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {

  ResetHyperparameters(ref_vector<Eigen::VectorXd>(inputs.begin()+1, inputs.end()));

  Eigen::VectorXd delta = inputs.at(0).get() - mean;

  return -0.5*varSize*std::log(2.0*M_PI) - 0.5*LogDeterminant() - 0.5 * delta.dot( ApplyPrecision(delta).col(0) );
}

Eigen::VectorXd GaussianBase::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) {

  ResetHyperparameters(ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.end()));

  Eigen::VectorXd z = RandomGenerator::GetNormal(mean.rows());

  return mean + ApplyCovSqrt(z);
}

unsigned int GaussianBase::Dimension() const {
  return mean.rows();
}

void GaussianBase::SetMean(Eigen::VectorXd const& newMu)
{
  assert(newMu.rows() == mean.rows());
  mean = newMu;
}

Eigen::VectorXd GaussianBase::GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  Eigen::VectorXd delta = inputs.at(0).get() - mean;
  if(wrt==0){
    return -1.0 * ApplyPrecision(delta);
  }else{
    std::cerr << "ERROR: Gradient wrt mean and covariance has not been implemented." << std::endl;
    assert(false);
  }
}
