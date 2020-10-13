#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Modeling;
using namespace muq::Utilities;

InverseGamma::InverseGamma(double       alphaIn,
                           double       betaIn) : Distribution(1),
                                                  alpha(alphaIn),
                                                  beta(betaIn),
                                                  logConst(alphaIn * std::log(betaIn) - std::lgamma(alphaIn))
                                                  {};

double InverseGamma::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  double x = inputs.at(0)(0);

  if(x<std::numeric_limits<double>::epsilon())
    return -1.0*std::numeric_limits<double>::infinity();

  return logConst + (-alpha-1.0)*std::log(x) - beta / x;
}


Eigen::VectorXd InverseGamma::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  Eigen::VectorXd output(1);
  output(0) = 1.0/RandomGenerator::GetGamma(alpha,1.0/beta);
  return output;
}
