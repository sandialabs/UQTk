#include "MUQ/Modeling/LinearAlgebra/KroneckerProductOperator.h"

#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

using namespace muq::Modeling;

KroneckerProductOperator::KroneckerProductOperator(std::shared_ptr<LinearOperator> Ain,
                                                   std::shared_ptr<LinearOperator> Bin) : LinearOperator(Ain->rows()*Bin->rows(), Ain->cols()*Bin->cols()), A(Ain), B(Bin)
{

}

Eigen::MatrixXd KroneckerProductOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    Eigen::MatrixXd output(nrows, x.cols());
    
    for(int i=0; i<x.cols(); ++i)
    {
        Eigen::VectorXd xVec = x.col(i);
        Eigen::Map<Eigen::MatrixXd> xMat(xVec.data(), B->cols(), A->cols());
        Eigen::Map<Eigen::MatrixXd> bMat(&output(0,i), B->rows(), A->rows());

        bMat = A->Apply( B->Apply( xMat ).transpose() ).transpose();
    }

    return output;
}



Eigen::MatrixXd KroneckerProductOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output(ncols, x.cols());
    
    for(int i=0; i<x.cols(); ++i)
    {
        Eigen::VectorXd xVec = x.col(i);
        Eigen::Map<const Eigen::MatrixXd> xMat(xVec.data(), B->rows(), A->rows());
        Eigen::Map<Eigen::MatrixXd> bMat(&output(0,i), B->cols(), A->cols());

        bMat = A->ApplyTranspose( B->ApplyTranspose( xMat ).transpose() ).transpose();
    }

    return output;
}


std::shared_ptr<LinearOperator> muq::Modeling::KroneckerSum(std::shared_ptr<LinearOperator> A,
                                                             std::shared_ptr<LinearOperator> B)
{
    auto part1 = std::make_shared<KroneckerProductOperator>(A, std::make_shared<IdentityOperator>(B->cols()));
    auto part2 = std::make_shared<KroneckerProductOperator>(std::make_shared<IdentityOperator>(A->rows()), B);

    return std::make_shared<SumOperator>(part1, part2);
}


Eigen::MatrixXd muq::Modeling::KroneckerProduct(Eigen::Ref<const Eigen::MatrixXd> const& A,
                                                 Eigen::Ref<const Eigen::MatrixXd> const& B)
{
    Eigen::MatrixXd output(A.rows()*B.rows(), A.cols()*B.cols());
    for(int j=0; j<A.cols(); ++j)
    {
        for(int i=0; i<A.rows(); ++i)
        {
            output.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
        }
    }
    return output;
}
