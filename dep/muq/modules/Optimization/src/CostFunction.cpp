#include "MUQ/Optimization/CostFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

CostFunction::CostFunction(Eigen::VectorXi const& inputSizes) :
  ModPiece(inputSizes, Eigen::VectorXi::Ones(1)) {}

CostFunction::~CostFunction() {}

void CostFunction::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  outputs.resize(1);
  outputs.at(0) = Eigen::VectorXd::Constant(1, CostImpl(input));
}

void CostFunction::GradientImpl(unsigned int const outputDimWrt,
                                unsigned int const inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd const& sensitivity) {
  GradientImpl(inputDimWrt, input, sensitivity);
}

void CostFunction::GradientImpl(unsigned int const inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd const& sensitivity) {
  ModPiece::GradientImpl(0, inputDimWrt, input, sensitivity);
}

double
CostFunction::Cost(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {
    return Evaluate(input).at(0) (0);
}

Eigen::VectorXd const&
CostFunction::Gradient(unsigned int const inputDimWrt,
                       std::vector<Eigen::VectorXd> const& input,
                       Eigen::VectorXd const& sensitivity) {
  return ModPiece::Gradient(0, inputDimWrt, input, sensitivity);
}


Eigen::MatrixXd
CostFunction::Hessian(unsigned int const inputDimWrt,
                      std::vector<Eigen::VectorXd> const& input) {
  return HessianByFD(inputDimWrt, input);
}


Eigen::MatrixXd
CostFunction::HessianByFD(unsigned int const inputDimWrt,
                          std::vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd sensitivity = Eigen::VectorXd::Ones(1);

  Eigen::VectorXd f0 = ModPiece::Gradient(0, inputDimWrt, input, sensitivity);
  Eigen::VectorXd f;

  double eps;
  Eigen::VectorXd newInput= input.at(inputDimWrt);
  std::vector<Eigen::VectorXd> newInputVec = input;

  Eigen::MatrixXd hes(inputSizes(inputDimWrt), inputSizes(inputDimWrt));
  for (int i=0; i<inputSizes(inputDimWrt); ++i) {
    eps = std::max(1.0e-8,
                   1.0e-10*std::abs(input.at(inputDimWrt)(i)));

    newInput(i) = input.at(inputDimWrt)(i) + eps;
    newInputVec.at(inputDimWrt) = std::cref(newInput);

    f = ModPiece::Gradient(0, inputDimWrt, newInputVec, sensitivity);

    hes.col(i) = (f-f0)/eps;

    newInput(i) = input.at(inputDimWrt)(i);

  }

  return hes;

}


Eigen::MatrixXd
CostFunction::ApplyHessian(unsigned int const inputDimWrt,
                           std::vector<Eigen::VectorXd> const& input,
                           Eigen::VectorXd const& vec) {

  return Hessian(inputDimWrt, input)*vec;

}
