#include "MUQ/Optimization/ModPieceCostFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

ModPieceCostFunction::ModPieceCostFunction(std::shared_ptr<ModPiece> cost) : CostFunction(cost->inputSizes), cost(cost) {
  // can only have one output of size one
  assert(cost->outputSizes.size()==1);
  assert(cost->outputSizes(0)==1);
}

double ModPieceCostFunction::CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {
  assert(cost);
  return cost->Evaluate(input).at(0) (0);
}

void ModPieceCostFunction::GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  assert(cost);
  gradient = cost->Gradient(0, inputDimWrt, input, sensitivity);
}
