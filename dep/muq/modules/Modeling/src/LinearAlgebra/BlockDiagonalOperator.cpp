#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"


using namespace muq::Modeling;


BlockDiagonalOperator::BlockDiagonalOperator(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn) : LinearOperator(SumRows(blocksIn), SumCols(blocksIn)), blocks(blocksIn){};


Eigen::MatrixXd BlockDiagonalOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    assert(x.rows() == ncols);
    
    int currRow = 0;
    int currCol = 0;

    Eigen::MatrixXd output(nrows,x.cols());
    
    for(int i=0; i<blocks.size(); ++i)
    {
        output.block(currRow, 0, blocks.at(i)->rows(), x.cols()) = blocks.at(i)->Apply( x.block(currCol, 0, blocks.at(i)->cols(), x.cols()) );
        currRow += blocks.at(i)->rows();
        currCol += blocks.at(i)->cols();
    }

    return output;
}

Eigen::MatrixXd BlockDiagonalOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{

    assert(x.rows() == nrows);
    
    int currRow = 0;
    int currCol = 0;

    Eigen::MatrixXd output(ncols,x.cols());
    
    for(int i=0; i<blocks.size(); ++i)
    {
        output.block(currRow, 0, blocks.at(i)->cols(), x.cols()) = blocks.at(i)->ApplyTranspose( x.block(currCol, 0, blocks.at(i)->rows(), x.cols()) );
        currRow += blocks.at(i)->cols();
        currCol += blocks.at(i)->rows();
    }

    return output;
}

Eigen::MatrixXd BlockDiagonalOperator::GetMatrix()
{

    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nrows,ncols);
    int currRow = 0;
    int currCol = 0;

    for(int i=0; i<blocks.size(); ++i)
    {
        output.block(currRow, currCol, blocks.at(i)->rows(), blocks.at(i)->cols()) = blocks.at(i)->GetMatrix();
        currRow += blocks.at(i)->rows();
        currCol += blocks.at(i)->cols();
    }

    return output;
}


int BlockDiagonalOperator::SumRows(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn)
{
    int sum = 0;
    for(auto& block : blocksIn)
        sum += block->rows();
    return sum;
}

int BlockDiagonalOperator::SumCols(std::vector<std::shared_ptr<LinearOperator>> const& blocksIn)
{
    int sum = 0;
    for(auto& block : blocksIn)
        sum += block->cols();
    return sum;
}
