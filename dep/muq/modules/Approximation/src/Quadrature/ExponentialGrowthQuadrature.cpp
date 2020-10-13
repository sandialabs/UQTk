#include "MUQ/Approximation/Quadrature/ExponentialGrowthQuadrature.h"

using namespace muq::Approximation;

ExponentialGrowthQuadrature::ExponentialGrowthQuadrature(std::shared_ptr<Quadrature> const& quadIn) : Quadrature(1),
                                                                                                      otherQuad(quadIn)
{
  assert(quadIn->Dim()==1);
}

void ExponentialGrowthQuadrature::Compute(unsigned int index)
{
  otherQuad->Compute(std::exp2(index));
}

unsigned int ExponentialGrowthQuadrature::Exactness(unsigned int index) const
{
  return otherQuad->Exactness(std::exp2(index));
}

Eigen::MatrixXd const& ExponentialGrowthQuadrature::Points() const
{
  return otherQuad->Points();
}

Eigen::VectorXd const& ExponentialGrowthQuadrature::Weights() const
{
  return otherQuad->Weights();
}
