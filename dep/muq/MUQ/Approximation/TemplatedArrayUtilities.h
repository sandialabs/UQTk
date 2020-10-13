
#ifndef TEMPLATEDARRAYUTILITIES_H_
#define TEMPLATEDARRAYUTILITIES_H_

#include <Eigen/Dense>
#include <vector>

namespace muq
{
namespace Approximation
{

template<typename MatType>
class MatrixBlock;

/** 

Templated function to get the shape of an array.  In Eigen parlance, GetShape(mat,0) is the same as mat.rows() and GetShape(mat,1) is the same as mat.cols().

This function is overloaded to handle general types of arrays.  Currently, Kokkos::View<double**> and Eigen matrices are supported.

*/
template<typename MatrixType>
unsigned GetShape(MatrixType const& mat, unsigned dim)
{
    return mat.dimension(dim);
}

template<typename ScalarType, int rows, int cols>
unsigned GetShape(Eigen::Matrix<ScalarType, rows, cols> const& mat, unsigned dim)
{
    assert(dim<2);
    return dim==0 ? mat.rows() : mat.cols();
}

template<typename Derived>
unsigned GetShape(MatrixBlock<Derived> const& mat, unsigned dim)
{
    assert(dim<2);
    return dim==0 ? mat.rows() : mat.cols();
}


template<typename Derived>
unsigned GetShape(Eigen::Ref<Derived> const& mat, unsigned dim)
{
    assert(dim<2);
    return dim==0 ? mat.rows() : mat.cols();
}

/** @brief Calculates the distance squared between two points defined by vectors v1 and v2. 
    @details Assumes the vectors are the same size and recursively compute the squared distance
             between them.  The recursion is used for numerical accuracy.
*/
template<typename VectorType1, typename VectorType2>
double CalcSquaredDistance(VectorType1 const& v1, VectorType2 const& v2, int startDim=0, int endDim=-1)
{
    if(endDim==-1)
    {
	endDim = GetShape(v1,0);
    }

    const int dim = endDim-startDim;
    const int minDim = 10;
    
    // If the dimension is small enough, just compute the some with a for loop
    if(dim<minDim)
    {
	double output = 0.0;
	for(int i=0; i<dim; ++i)
	{
	    output += std::pow(v1(startDim+i)-v2(startDim+i), 2.0);
	}
	return output;
    }
    else
    {

	int midDim = startDim + std::floor(0.5*dim);
	return CalcSquaredDistance(v1,v2, startDim, midDim) + CalcSquaredDistance(v1,v2, midDim, endDim);
    }

}


/** Calculates the distance between two points defined by vectors v1 and v2. */
template<typename VectorType1, typename VectorType2>
double CalcDistance(VectorType1 const& v1, VectorType2 const& v2)
{
    // Make sure the vectors are the same size
    const int dim = GetShape(v1,0);
    assert(dim==GetShape(v2,0));

    return std::sqrt(CalcSquaredDistance(v1,v2));
}

/** 

This class is used like mat.col(colNum) in Eigen, but can be used with more general matrix types (e.g., Kokkos::View<double**>).  Based on a two dimensional matrix (anything implementing operator()(int,int), this class provides a mechanism for accessing elements in a particular column.

This class is usually not constructed directly, but is created using the "ColumnSlice" function, allowing the compiler to detect the correct template type.

 */
template<typename MatType>
class ColumnSlice
{
public:
    //ColumnSlice(MatType const& matrixIn, unsigned colIn) : col(colIn), matrix(matrixIn){};
    ColumnSlice(MatType& matrixIn, unsigned colIn) : col(colIn), matrix(matrixIn){};

    double operator()(unsigned row) const{return matrix(row, col);};
    double operator()(unsigned row, unsigned col2) const{assert(col2==0); return matrix(row,col);};

    double& operator()(unsigned row){return matrix(row, col);};
    double& operator()(unsigned row, unsigned col2){assert(col2==0); return matrix(row,col);};
    
    unsigned dimension(unsigned dim) const
    {
	if(dim>0)
	    return 1;
	else
	    return GetShape(matrix,0);
    };

    unsigned rows() const
    {
        return GetShape(matrix,0);
    }

    unsigned cols() const
    {
        return 1;
    }

    operator Eigen::VectorXd() const
    {
      return eval();
    };
    
    Eigen::VectorXd eval() const
    {
        Eigen::VectorXd output(rows());
        for(int i=0; i<rows(); ++i)
            output(i) = (*this)(i);
        return output;
    };
private:
    const unsigned   col;
    MatType& matrix;
};

template<typename VecType>
class VectorSlice
{
public:
    VectorSlice(VecType const& vectorIn, std::vector<unsigned> const& indsIn) : vector(vectorIn), inds(indsIn){};

    double  operator()(unsigned row) const{return vector(inds.at(row));};
    double& operator()(unsigned row){return vector(inds.at(row));};

    unsigned dimension(unsigned dim) const
    {
	if(dim>0)
	    return 1;
	else
	    return inds.size();
    }
    
private:
    VecType const& vector;
    std::vector<unsigned> const& inds;
};

template<typename MatType>
class MatrixBlock
{
public:
    MatrixBlock(MatType & matrixIn,
		unsigned startRowIn,
		unsigned startColIn,
		unsigned numRowsIn,
		unsigned numColsIn) : startRow(startRowIn), startCol(startColIn), numRows(numRowsIn), numCols(numColsIn), matrix(matrixIn)
    {
	assert(GetShape(matrix,0)>=startRow+numRows);
	assert(GetShape(matrix,1)>=startCol+numCols);
    };

    double operator() (unsigned row, unsigned col) const{return matrix(startRow + row, startCol + col);};
    double &operator()(unsigned row, unsigned col){return matrix(startRow + row, startCol + col);};

    template<typename Derived>
    MatrixBlock& operator=(Eigen::DenseBase<Derived> const& otherMat)
    {
	assert(otherMat.rows()==numRows);
	assert(otherMat.cols()==numCols);

	for(unsigned j=0; j<numCols; ++j)
	{
	    for(unsigned i=0; i<numRows; ++i)
	    {
		matrix(startRow + i, startCol + j) = otherMat(i,j);
	    }
	}
	return *this;
    }
    
    unsigned rows() const{return numRows;};
    unsigned cols() const{return numCols;};
    
private:
    
    const unsigned startRow, startCol, numRows, numCols;
    MatType& matrix;
};

/** 

Grab a particular column of a matrix and return as an instance of the "ColumnSlice" class. 

*/
template<typename MatType>
ColumnSlice<MatType> GetColumn(MatType& matrix, unsigned col)
{
    return ColumnSlice<MatType>(matrix, col);
}

template<typename MatType>
VectorSlice<MatType> GetSlice(MatType& matrix, std::vector<unsigned> const& inds)
{
    return VectorSlice<MatType>(matrix, inds);
}


template<typename MatType>
MatrixBlock<MatType> GetBlock(MatType & matrix, unsigned rowStart, unsigned colStart, unsigned numRows, unsigned numCols)
{
    return MatrixBlock<MatType>(matrix, rowStart, colStart, numRows, numCols);
}


} // namespace Approximation

} // namespace muq

#endif 
