#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

ConstantVector::ConstantVector(std::vector<Eigen::VectorXd> const& outs) : ModPiece(Eigen::VectorXi(), OutSizes(outs))
{
  outputs.resize(outs.size());
  for( unsigned int i=0; i<outs.size(); ++i ) {
    outputs.at(i) = outs[i];
  }
}


ConstantVector::ConstantVector(Eigen::VectorXd const& valIn) : ModPiece(Eigen::VectorXi(), valIn.size()*Eigen::VectorXi::Ones(1))
{
  outputs.resize(1);
  outputs.at(0) = valIn;
}

Eigen::VectorXi ConstantVector::OutSizes(std::vector<Eigen::VectorXd> const& outs) {
  Eigen::VectorXi oSizes(outs.size());
  for( unsigned int i=0; i<outs.size(); ++i ) {
    oSizes(i) = outs[i].size();
  }

  return oSizes;
}

void ConstantVector::SetValue(Eigen::VectorXd const& valIn)
{
  if(valIn.size() != outputSizes(0))
    throw muq::WrongSizeError("In ConstantVector::SetValue, new vector has size " + std::to_string(valIn.size()) + ", but expected a size of " + std::to_string(outputSizes(0)) + ".");

  outputs.at(0) = valIn;
}

void ConstantVector::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs){
  return;
}
