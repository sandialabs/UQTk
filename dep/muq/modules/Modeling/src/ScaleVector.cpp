#include "MUQ/Modeling/ScaleVector.h"

using namespace muq::Modeling;

ScaleVector::ScaleVector(double const scale, unsigned int const dim) : ModPiece(Eigen::VectorXi::Constant(1, dim), Eigen::VectorXi::Constant(1, dim)), scale(scale) {}

void ScaleVector::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  outputs.resize(1);
  outputs[0] = scale*inputs[0];
}

void ScaleVector::JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) {
  jacobian = scale*Eigen::MatrixXd::Identity(outputSizes(0), inputSizes(0));
}

void ScaleVector::GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) {
  gradient = scale*sens;
}

void ScaleVector::ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) {
  jacobianAction = scale*targ;
}
