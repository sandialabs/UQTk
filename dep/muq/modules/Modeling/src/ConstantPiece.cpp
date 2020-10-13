#include "MUQ/Modeling/ConstantPiece.h"

using namespace muq::Modeling;

ConstantPiece::ConstantPiece(std::vector<boost::any> const& outs) : WorkPiece(0, Types(outs)) {
  // the outputs will not change inside evaluate so we should not clear them
  clearOutputs = false;

  // populate the vector of outputs
  outputs.resize(numOutputs);
  std::copy(outs.begin(), outs.end(), outputs.begin());
}

ConstantPiece::ConstantPiece() : WorkPiece(0, -1) {
  // the outputs will not change inside evaluate so we should not clear them
  clearOutputs = false;
}

void ConstantPiece::SetOutputs(std::vector<boost::any> const& outs) {
  // clear the output types
  outputTypes.clear();

  // set the output types
  outputTypes = Types(Types(outs));

  // populate the vector of outputs
  outputs.resize(outs.size());
  std::copy(outs.begin(), outs.end(), outputs.begin());
}

void ConstantPiece::SetOutputs() {
  // clear the output types
  outputTypes.clear();

  // clear the vector of outputs
  outputs.clear();
}

// the outputs are already set and not cleared so don't do anything
void ConstantPiece::EvaluateImpl(ref_vector<boost::any> const& inputs) {}
