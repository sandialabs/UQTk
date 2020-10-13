#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

std::shared_ptr<Density> Distribution::AsDensity()
{
  return std::make_shared<Density>(shared_from_this());
}

std::shared_ptr<RandomVariable> Distribution::AsVariable()
{
  return std::make_shared<RandomVariable>(shared_from_this());
}

ref_vector<const Eigen::VectorXd> Distribution::ToRefVector(std::vector<Eigen::VectorXd> const& vec) const {

  ref_vector<const Eigen::VectorXd> refs;
  refs.reserve(vec.size());

  // populate the input vector
  for(int i=0; i<vec.size(); ++i)
    refs.push_back(std::cref(vec.at(i)));

  return refs;
}

double Distribution::LogDensity(ref_vector<Eigen::VectorXd> const& inputs) {
  return LogDensityImpl(inputs);
}


Eigen::VectorXd Distribution::Sample(ref_vector<Eigen::VectorXd> const& inputs) {
  return SampleImpl(inputs);
}

Eigen::VectorXd Distribution::Sample() {
  return Sample(ref_vector<Eigen::VectorXd>());
}

Eigen::VectorXd Distribution::GradLogDensity(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  assert(wrt<inputs.size());
  return GradLogDensityImpl(wrt, inputs);
}

Eigen::VectorXd Distribution::GradLogDensityImpl(unsigned int                       wrt,
                                                 ref_vector<Eigen::VectorXd> const& inputs)
{
  // Default to finite difference
  ref_vector<Eigen::VectorXd> newInputs = inputs;
  Eigen::VectorXd newIn = newInputs.at(wrt).get();
  newInputs.at(wrt) = std::cref(newIn);

  const double f0 = LogDensity(newInputs);

  const int dim = inputs.at(wrt).get().size();

  Eigen::VectorXd output(dim);
  for(int i=0; i<dim; ++i){

    double eps = std::max(1e-8, 1e-10*std::abs(inputs.at(wrt).get()(i)));
    newIn(i) += eps;

    double newF = LogDensity(newInputs);
    output(i) = (newF-f0)/eps;
    newIn(i) = inputs.at(wrt).get()(i);
  }

  return output;
}
