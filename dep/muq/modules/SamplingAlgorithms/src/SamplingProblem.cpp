#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingProblem::SamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& targetIn) : AbstractSamplingProblem(targetIn->inputSizes),
                                                                                      target(targetIn) {}



double SamplingProblem::LogDensity(unsigned int const t, std::shared_ptr<SamplingState> const& state, AbstractSamplingProblem::SampleType type) {
  assert(target);

  return target->Evaluate(state->state).at(0)(0);
}

Eigen::VectorXd SamplingProblem::GradLogDensity(std::shared_ptr<SamplingState> const& state,
                                                unsigned                       blockWrt)
{
  return target->Gradient(0,blockWrt, state->state, Eigen::VectorXd::Ones(1).eval());
}

std::vector<int> SamplingProblem::GetBlockSizes(std::shared_ptr<ModPiece> const& target)
{
  int numBlocks = target->inputSizes.size();

  std::vector<int> output(numBlocks);
  for(int i=0; i<numBlocks; ++i)
    output.at(i) = target->inputSizes(i);

  return output;
}
