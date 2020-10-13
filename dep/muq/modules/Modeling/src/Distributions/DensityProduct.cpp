#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Utilities/AnyHelpers.h"


using namespace muq::Modeling;
using namespace muq::Utilities;

DensityProduct::DensityProduct(int numPiecesIn) : DensityBase(Eigen::VectorXi::Ones(numPiecesIn))
{
  numInputs = numPiecesIn;
}


double DensityProduct::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  double sum = 0.0;
  for(int i=0; i<inputs.size(); ++i)
    sum += inputs.at(i).get()(0);

  return sum;
}

Eigen::VectorXd DensityProduct::GradLogDensityImpl(unsigned int wrt,
                                                   ref_vector<Eigen::VectorXd> const& inputs)
{
  return Eigen::VectorXd::Ones(1,1);
}
