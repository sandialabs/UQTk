#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"


using namespace muq::Modeling;


SumOperator::SumOperator(std::shared_ptr<LinearOperator> Ain,
                         std::shared_ptr<LinearOperator> Bin) : LinearOperator(Ain->rows(), Ain->cols()), A(Ain), B(Bin)
{
    assert(A->rows()==B->rows());
    assert(A->cols()==B->cols());
}


Eigen::MatrixXd SumOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    return A->Apply(x) + B->Apply(x);
}


Eigen::MatrixXd SumOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    return A->ApplyTranspose(x) + B->ApplyTranspose(x);
}
