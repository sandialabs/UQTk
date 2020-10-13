#include "MUQ/Modeling/CombineVectors.h"

using namespace muq::Modeling;

CombineVectors::CombineVectors(Eigen::VectorXi const& inputSizes) : ModPiece(inputSizes, Eigen::VectorXi::Constant(1, inputSizes.sum())) {}

void CombineVectors::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  Eigen::VectorXd out(outputSizes(0));
  unsigned int ind = 0;
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    out.segment(ind, inputs[i].get().size()) = inputs[i].get();
    ind += inputs[i].get().size();
  }

  outputs.resize(1);
  outputs[0] = out;
}

void CombineVectors::JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) {
  assert(outwrt==0);

  jacobian = Eigen::MatrixXd::Zero(outputSizes(0), inputSizes(inwrt));
  const unsigned int ind = inputSizes.segment(0, inwrt).sum();
  jacobian.block(ind, 0, inputSizes(inwrt), inputSizes(inwrt)) = Eigen::MatrixXd::Identity(inputSizes(inwrt), inputSizes(inwrt));
}

void CombineVectors::GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) {
  assert(outwrt==0);

  const unsigned int ind = inputSizes.segment(0, inwrt).sum();
  gradient = sens.segment(ind, inputSizes(inwrt));
}

void CombineVectors::ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) {
  assert(outwrt==0);

  const unsigned int ind = inputSizes.segment(0, inwrt).sum();
  jacobianAction = Eigen::VectorXd::Zero(outputSizes(0));
  jacobianAction.segment(ind, inputSizes(inwrt)) = targ;
}
