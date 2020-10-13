#include "MUQ/Modeling/Distributions/RandomVariable.h"

using namespace muq::Modeling;


RandomVariable::RandomVariable(std::shared_ptr<Distribution> distIn) : Distribution(distIn->varSize, distIn->hyperSizes),
                                                                       ModPiece(distIn->hyperSizes,
                                                                                distIn->varSize*Eigen::VectorXi::Ones(1)),
                                                                       dist(distIn)
{
  assert(dist);
}


void RandomVariable::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  outputs.resize(1);
  outputs.at(0) = dist->Sample(inputs);
}

double RandomVariable::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->LogDensity(inputs);
};

Eigen::VectorXd RandomVariable::GradLogDensityImpl(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->GradLogDensity(wrt, inputs);
}

Eigen::VectorXd RandomVariable::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->Sample(inputs);
}
