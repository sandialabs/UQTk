#include "MUQ/Modeling/SumPiece.h"

using namespace muq::Modeling;



SumPiece::SumPiece(unsigned int dim, unsigned int numInputs) : ModPiece(dim*Eigen::VectorXi::Ones(numInputs),
                                                                        dim*Eigen::VectorXi::Ones(1))
{
  assert(numInputs>1);
}



void SumPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input)
{
  outputs.resize(1);
  outputs.at(0) = input.at(0).get();

  for(unsigned int i=1; i<input.size(); ++i)
    outputs.at(0) += input.at(i).get();
}


void SumPiece::GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity)
{
  gradient = sensitivity;
}

void SumPiece::JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = Eigen::MatrixXd::Identity(inputSizes(inputDimWrt),inputSizes(inputDimWrt));
}

void SumPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec)
{
  jacobianAction = vec;
}

// void SumPiece::ApplyHessianImpl(unsigned int                const  outWrt,
//                                    unsigned int                const  inWrt1,
//                                    unsigned int                const  inWrt2,
//                                    ref_vector<Eigen::VectorXd> const& input,
//                                    Eigen::VectorXd             const& sens,
//                                    Eigen::VectorXd             const& vec)
// {
//   hessAction = Eigen::VectorXd::Zero(vec.size());
// }
