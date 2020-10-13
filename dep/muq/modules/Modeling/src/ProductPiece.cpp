#include "MUQ/Modeling/ProductPiece.h"

using namespace muq::Modeling;

ProductPiece::ProductPiece(unsigned int vectorDim) : ModPiece(Eigen::VectorXi::Constant(1,vectorDim), Eigen::VectorXi::Ones(1)){};

void ProductPiece::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs)
{
  outputs.resize(0);
  outputs.at(0).resize(1);
  outputs.at(0)(0) = inputs.at(0).get().prod();
}
