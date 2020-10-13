#include "MUQ/Modeling/LinearAlgebra/ConcatenateOperator.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;

ConcatenateOperator::ConcatenateOperator(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                                         const int                                           rowColIn) : LinearOperator(GetRows(opsIn, rowColIn), GetCols(opsIn, rowColIn)), ops(opsIn), rowCol(rowColIn)
{
    CheckSizes();
};


/** Apply the linear operator to a vector */
Eigen::MatrixXd ConcatenateOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nrows, x.cols());

    if(rowCol==0){
        int currRow = 0;
        for(int i=0; i<ops.size(); ++i){
            output.block(currRow,0,ops.at(i)->rows(), x.cols()) = ops.at(i)->Apply(x);
            currRow += ops.at(i)->rows();
        }
    }else{
        int currRow = 0;
        for(int i=0; i<ops.size(); ++i){
            output += ops.at(i)->Apply( x.block(currRow, 0, ops.at(i)->cols(), x.cols()) );
            currRow += ops.at(i)->cols();
        }
    }

    return output;
}


/** Apply the transpose of the linear operator to a vector. */
Eigen::MatrixXd ConcatenateOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(ncols, x.cols());

    if(rowCol==0){
        int currRow = 0;
        for(int i=0; i<ops.size(); ++i){
            output += ops.at(i)->ApplyTranspose( x.block(currRow, 0, ops.at(i)->rows(), x.cols()) );
            currRow += ops.at(i)->rows();
        }
    }else{
        int currRow = 0;
        for(int i=0; i<ops.size(); ++i){
            output.block(currRow,0, ops.at(i)->cols(), x.cols()) = ops.at(i)->ApplyTranspose(x);
            currRow += ops.at(i)->cols();
        }
    }

    return output;

}


 std::shared_ptr<ConcatenateOperator> ConcatenateOperator::VStack(std::shared_ptr<LinearOperator> Ain,
                                                                  std::shared_ptr<LinearOperator> Bin)
 {
     std::vector<std::shared_ptr<LinearOperator>> temp{Ain, Bin};
     return std::make_shared<ConcatenateOperator>(temp, 0);
 }

std::shared_ptr<ConcatenateOperator> ConcatenateOperator::HStack(std::shared_ptr<LinearOperator> Ain,
                                                                 std::shared_ptr<LinearOperator> Bin)
{
    std::vector<std::shared_ptr<LinearOperator>> temp{Ain, Bin};
    return std::make_shared<ConcatenateOperator>(temp,1);
}


Eigen::MatrixXd ConcatenateOperator::GetMatrix()
{

    Eigen::MatrixXd output(nrows, ncols);
    if(rowCol==0){
        int currRow = 0;
        for(int i=0; i<ops.size(); ++i){
            output.block(currRow, 0, ops.at(i)->rows(), ncols) = ops.at(i)->GetMatrix();
            currRow += ops.at(i)->rows();
        }
    }else{
        int currCol = 0;
        for(int i=0; i<ops.size(); ++i){
            output.block(0,currCol, nrows, ops.at(i)->cols()) = ops.at(i)->GetMatrix();
            currCol += ops.at(i)->cols();
        }
    }
    return output;
}

int ConcatenateOperator::GetRows(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                                 const int                                           rowColIn)
{
    assert(opsIn.size()>0);

    if(rowColIn==0){
        int count = 0;
        for(auto& op : opsIn)
            count += op->rows();
        return count;
    }else{
        return opsIn.at(0)->rows();
    }
}


int ConcatenateOperator::GetCols(std::vector<std::shared_ptr<LinearOperator>> const& opsIn,
                                 const int                                           rowColIn)
{
    assert(opsIn.size()>0);

    if(rowColIn==0){
        return opsIn.at(0)->cols();
    }else{
        int count = 0;
        for(auto& op : opsIn)
            count += op->cols();
        return count;
    }
}


void ConcatenateOperator::CheckSizes()
{   
    if(rowCol==0){
        const int fixedCols = ops.at(0)->cols();
        for(int i=1; i<ops.size(); ++i){
            if(fixedCols != ops.at(i)->cols())
                throw muq::WrongSizeError("In ConcatenateOperator: Cannot vertically stack operators with different number of columns.  Matrix A has " + std::to_string(fixedCols) + " columns but matrix B has " + std::to_string(ops.at(i)->cols()) + " columns.");
        }

    }else{
        const int fixedRows = ops.at(0)->rows();
        for(int i=1; i<ops.size(); ++i){
            if(fixedRows != ops.at(i)->rows())
                throw muq::WrongSizeError("In ConcatenateOperator: Cannot horizontally stack operators with different number of rows.  Matrix A has " + std::to_string(fixedRows) + " rows but matrix B has " + std::to_string(ops.at(i)->rows()) + " rows.");
        }
    }
}
