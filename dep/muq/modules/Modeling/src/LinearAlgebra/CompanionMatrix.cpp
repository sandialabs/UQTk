#include "MUQ/Modeling/LinearAlgebra/CompanionMatrix.h"

using namespace muq::Modeling;



Eigen::MatrixXd CompanionMatrix::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    assert(x.rows() == ncols);
    
    Eigen::MatrixXd output(nrows, x.cols());
    
    output.block(0, 0, nrows-1, x.cols()) = x.block(1,0,nrows-1, x.cols());
    output.row(nrows-1) = lastRow.transpose()*x;

    return output;
}


Eigen::MatrixXd CompanionMatrix::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    assert(x.rows() == nrows);
    
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(ncols, x.cols());

    output.block(1, 0, ncols-1, x.cols()) = x.block(0,0,ncols-1,x.cols());
    output += lastRow * x.row(x.rows()-1);
    
    return output;
}


Eigen::MatrixXd CompanionMatrix::GetMatrix()
{

    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nrows, ncols);
    output.block(0,1,nrows-1,nrows-1) = Eigen::MatrixXd::Identity(nrows-1,nrows-1);
    output.row(nrows-1) = lastRow;

    return output;
}
