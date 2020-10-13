#include "MUQ/Modeling/SplitVector.h"

using namespace muq::Modeling;

SplitVector::SplitVector(Eigen::VectorXi const& ind, Eigen::VectorXi const& size, unsigned int const insize) : ModPiece(Eigen::VectorXi::Constant(1, insize), size), ind(ind), size(size) {
  assert(ind.size()==size.size());
  assert(size.sum()<=insize);
  assert(ind.maxCoeff()<insize);
}

void SplitVector::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  const Eigen::VectorXd& in = inputs[0];

  outputs.resize(ind.size());
  for( unsigned int i=0; i<ind.size(); ++i ) {
    outputs[i] = in.segment(ind(i), size(i)).eval();
  }
}

void SplitVector::JacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs) {
  assert(inwrt==0);
  jacobian = Eigen::MatrixXd::Zero(size(outwrt), inputSizes(0));
  jacobian.block(0, ind(outwrt), size(outwrt), size(outwrt)) += Eigen::MatrixXd::Identity(size(outwrt), size(outwrt)).eval();
}

void SplitVector::GradientImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& sens) {
  assert(inwrt==0);
  assert(sens.size()==size(outwrt));
  gradient = Eigen::VectorXd::Zero(inputSizes(0));
  gradient.segment(ind(outwrt), size(outwrt)) += sens;
}

void SplitVector::ApplyJacobianImpl(unsigned int const outwrt, unsigned int const inwrt, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& targ) {
  assert(inwrt==0);
  assert(targ.size()==inputSizes(0));
  jacobianAction = targ.segment(ind(outwrt), size(outwrt));
}
