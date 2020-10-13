#include "MUQ/Modeling/LinearAlgebra/ProductOperator.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

ProductOperator::ProductOperator(std::shared_ptr<LinearOperator> Ain,
                                 std::shared_ptr<LinearOperator> Bin) : LinearOperator(Ain->rows(), Bin->cols()), A(Ain), B(Bin)
{

    if(A->cols() != B->rows())
        throw muq::WrongSizeError("In ProductOperator: The number of columns in A (" + std::to_string(A->cols()) + ") must match the number of rows in B (" + std::to_string(B->rows()) + ")");
}

/** Apply the linear operator to a vector */
Eigen::MatrixXd ProductOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    return A->Apply(B->Apply(x));
}

/** Apply the transpose of the linear operator to a vector. */
Eigen::MatrixXd ProductOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    return B->ApplyTranspose(A->ApplyTranspose(x));
}

Eigen::MatrixXd ProductOperator::GetMatrix()
{
    return A->GetMatrix()*B->GetMatrix();
}
