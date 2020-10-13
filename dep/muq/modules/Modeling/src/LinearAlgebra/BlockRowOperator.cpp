#include "MUQ/Modeling/LinearAlgebra/BlockRowOperator.h"


using namespace muq::Modeling;


BlockRowOperator::BlockRowOperator(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn) : LinearOperator(blocksIn.at(0)->rows(), SumCols(blocksIn)), blocks(blocksIn)
{

    for(int i=0; i<blocks.size(); ++i)
        assert(blocks.at(i)->rows()==nrows);

};


Eigen::MatrixXd BlockRowOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    assert(x.rows() == ncols);
    
    int currCol = 0;

    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nrows,x.cols());
    
    for(int i=0; i<blocks.size(); ++i)
    {
        output += blocks.at(i)->Apply( x.block(currCol, 0, blocks.at(i)->cols(), x.cols()) );
        currCol += blocks.at(i)->cols();
    }

    return output;
}

Eigen::MatrixXd BlockRowOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    assert(x.rows() == nrows);
    
    int currRow = 0;
    int currCol = 0;

    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(ncols,x.cols());
    
    for(int i=0; i<blocks.size(); ++i)
    {
        output.block(currRow, 0, blocks.at(i)->cols(), x.cols()) = blocks.at(i)->ApplyTranspose( x );
        currRow += blocks.at(i)->cols();
        currCol += blocks.at(i)->rows();
    }

    return output;
}

Eigen::MatrixXd BlockRowOperator::GetMatrix()
{

    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nrows,ncols);
    int currCol = 0;

    for(int i=0; i<blocks.size(); ++i)
    {
        output.block(0, currCol, nrows, blocks.at(i)->cols()) = blocks.at(i)->GetMatrix();
        currCol += blocks.at(i)->cols();
    }

    return output;
}

int BlockRowOperator::SumCols(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn)
{
    int sum = 0;
    for(auto& block : blocksIn)
        sum += block->cols();
    return sum;
}
