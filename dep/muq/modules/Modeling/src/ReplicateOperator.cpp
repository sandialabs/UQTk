#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;


ReplicateOperator::ReplicateOperator(unsigned int vectorDim,
                                     unsigned int numRepeat) : ModPiece(int(vectorDim)*Eigen::VectorXi::Ones(1),
                                                                        int(vectorDim*numRepeat)*Eigen::VectorXi::Ones(1)),
                                                               numRepl(numRepeat)
{}

void ReplicateOperator::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs){
  outputs.resize(1);
  outputs.at(0).resize(outputSizes(0));

  unsigned int vectorDim = inputs.at(0).get().size();

  for(int i=0; i<numRepl; ++i)
    outputs.at(0).segment(i*vectorDim, vectorDim) = inputs.at(0).get();
}
