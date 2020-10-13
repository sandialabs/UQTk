#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"


using namespace muq::Modeling;

LinearOperator::LinearOperator(int rowsIn, int colsIn, int numInputCols) : muq::Modeling::ModPiece(colsIn*numInputCols*Eigen::VectorXi::Ones(1),
                                                                                                   rowsIn*numInputCols*Eigen::VectorXi::Ones(1)),
                                                         ncols(colsIn),
                                                         nrows(rowsIn)
{
};

void LinearOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = Apply(x);
};


void LinearOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y)
{
    assert(y.cols()==x.cols());
    y = ApplyTranspose(x);
};

Eigen::MatrixXd LinearOperator::GetMatrix()
{

    Eigen::MatrixXd output(nrows, ncols);
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Identity(ncols, ncols);

    for(int i=0; i<ncols; ++i)
        output.col(i) = Apply(rhs.col(i));

    return output;
}

void LinearOperator::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  outputs.resize(1);
  outputs.at(0) = Apply(input.at(0).get()).col(0);
}

void LinearOperator::GradientImpl(unsigned int                const  outputDimWrt,
                          unsigned int                const  inputDimWrt,
                          muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                          Eigen::VectorXd             const& sensitivity)
{
  gradient = ApplyTranspose(sensitivity);
}

void LinearOperator::JacobianImpl(unsigned int                const  outputDimWrt,
                          unsigned int                const  inputDimWrt,
                          muq::Modeling::ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = GetMatrix();
}

void LinearOperator::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                               unsigned int                const  inputDimWrt,
                               muq::Modeling::ref_vector<Eigen::VectorXd> const& input,
                               Eigen::VectorXd             const& vec)
{
  jacobianAction = Apply(vec);
}
