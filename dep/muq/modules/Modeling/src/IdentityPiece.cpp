#include "MUQ/Modeling/IdentityPiece.h"

using namespace muq::Modeling;

// Create a muq::Modeling::IdentityPiece with no fixed number of inputs and outputs and variable input/output types.
IdentityPiece::IdentityPiece() : WorkPiece() {}

// Create a muq::Modeling::IdentityPiece with fixed number of inputs/outputs and variable input/output types.
IdentityPiece::IdentityPiece(int const num) : WorkPiece(num, num) {}

// Create a muq::Modeling::IdentityPiece with a fixed number of inputs/outputs with specified types 
IdentityPiece::IdentityPiece(std::vector<std::string> const& types) : WorkPiece(types, types) {}

// Create a muq::Modeling::IdentityPiece where some of the inputs/outputs have specified types 
IdentityPiece::IdentityPiece(std::map<unsigned int, std::string> const& types) : WorkPiece(types, types) {}

// Create a muq::Modeling::IdentityPiece where some of the inputs/outputs have specified types the number of inputs/outputs is fixed
IdentityPiece::IdentityPiece(std::map<unsigned int, std::string> const& types, unsigned int const num) : WorkPiece(types, num, types, num) {}

void IdentityPiece::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // the number of outputs is the name as the number of inputs
  outputs.resize(inputs.size());

  // copy the value of the inputs
  for( unsigned int i=0; i<outputs.size(); ++i ) {
    outputs[i] = inputs.at(i).get();
  }
}
