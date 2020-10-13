#include "MUQ/Modeling/LinearAlgebra/AffineOperator.h"


using namespace muq::Modeling;

AffineOperator::AffineOperator(std::shared_ptr<LinearOperator> const& Ain,
                               Eigen::VectorXd                 const& bIn) : ModPiece(Ain->inputSizes,
                                                                                      Ain->outputSizes),
                                                                             A(Ain),
                                                                             b(bIn)
{
  assert(A->rows()==b.rows());
};

void AffineOperator::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  outputs.resize(1);
  outputs.at(0) = b + A->Apply(input.at(0).get()).col(0);
}

void AffineOperator::GradientImpl(unsigned int                const  outputDimWrt,
                                  unsigned int                const  inputDimWrt,
                                  muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                                  Eigen::VectorXd             const& sensitivity)
{
  gradient = A->ApplyTranspose(sensitivity);
}

void AffineOperator::JacobianImpl(unsigned int                const  outputDimWrt,
                                  unsigned int                const  inputDimWrt,
                                  muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = A->GetMatrix();
}

void AffineOperator::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                       unsigned int                const  inputDimWrt,
                                       muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& vec)
{
  jacobianAction = A->Apply(vec);
}
