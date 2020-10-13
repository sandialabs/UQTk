#ifndef LYAPUNOVSOLVER_H
#define LYAPUNOVSOLVER_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace muq
{
namespace Modeling
{

    /** @class LyaponovSolver
        @brief Solves the Lyapunov equation \f$A^TX + XA + Q = 0\f$ given square matrices \f$A\f$ and \f$Q\f$.  
        @details This class solves the continuous time Lyapunov equation given by \f[A^TX + XA + Q = 0,\f] where \f$A\f$ and \f$Q\f$ are square matrices.
                 In general, this class follows a similar form to the eigenvalue solvers in Eigen; the "compute" function is called with \f$A\f$ and \f$Q\f$ to compute the solution \f$X\f$, and then the solution \f$X\f$ is accessed using the "matrixX" function.   For example,

\code{.cpp}
Eigen::MatrixXd A;
// Fill in A here

Eigen::MatrixXd Q;
// Fill in Q here

LyapunovSolver<double> solver;
solver.compute(A,Q);

Eigen::MatrixXcd const& X = solver.matrixX();
\endcode

     */
    template<class ScalarType, int FixedRows=Eigen::Dynamic, int FixedCols=Eigen::Dynamic>
    class LyapunovSolver
    {
    public:

        typedef Eigen::Matrix<ScalarType, FixedRows, FixedCols> MatrixType;
        typedef Eigen::Matrix<std::complex<ScalarType>, FixedRows, FixedCols> ComplexMatrixType;
        
        LyapunovSolver& compute(MatrixType const& A, MatrixType const& C)
        {
            const int dim = A.rows();
            assert(A.rows()==A.cols());
            assert(C.rows()==C.cols());
            assert(A.rows()==C.rows());
            
            Eigen::ComplexSchur<MatrixType> schur;
            schur.compute(A);

            auto& Q = schur.matrixU();
            ComplexMatrixType S = schur.matrixT();

            ComplexMatrixType ctilde = Q.adjoint() * C.template cast<std::complex< ScalarType >>() * Q;

            X.resize(dim,dim);
            ComputeFromSchur(S, ctilde, X);

            X = (Q*X*Q.adjoint()).eval();

            return *this;
        };

        ComplexMatrixType const& matrixX() const{return X;};

    private:

        void ComputeFromSchur(Eigen::Ref<const ComplexMatrixType> const& S,
                              Eigen::Ref<ComplexMatrixType>              ctilde,
                              Eigen::Ref<ComplexMatrixType>              X)
        {
            const int size = X.rows();
            
            X(0,0) = -ctilde(0,0)/(S(0,0)+ std::conj(S(0,0)));

            if(size==1)
                return;
            
            ComplexMatrixType tempX = -ctilde.block(1,0,size-1,1) - X(0,0)*S.block(0,1,1,size-1).adjoint();
            (S.block(1,1,size-1,size-1).adjoint() + S(0,0)*ComplexMatrixType::Identity(size-1,size-1)).template triangularView<Eigen::Lower>().solveInPlace(tempX);

            X.block(1,0,size-1,1) = tempX;
            X.block(0,1,1,size-1) = X.block(1,0,size-1,1).adjoint();

            // Recursively call this function on the next block
            ctilde.block(1,1,size-1,size-1) += S.block(0,1,1,size-1).adjoint()*X.block(0,1,1,size-1) + X.block(1,0,size-1,1)*S.block(0,1,1,size-1);
                    
            ComputeFromSchur(S.block(1,1,size-1,size-1), ctilde.block(1,1,size-1,size-1) , X.block(1,1,size-1,size-1));
            
        }
    
        ComplexMatrixType X;
        
    };
}
}


#endif
