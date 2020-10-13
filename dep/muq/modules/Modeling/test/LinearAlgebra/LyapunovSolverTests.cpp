#include "MUQ/Modeling/LinearAlgebra/LyapunovSolver.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

using namespace muq::Modeling;


TEST(Utilities_LyapunovSolver, Diagonal)
{
    const int dim = 10;
    Eigen::MatrixXd A = -Eigen::MatrixXd::Identity(dim,dim);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(dim,dim);


    LyapunovSolver<double> solver;
    solver.compute(A,Q);

    auto X = solver.matrixX().real();

    for(int i=0; i<dim; ++i)
    {
        EXPECT_DOUBLE_EQ(0.5, X(i,i) );

        for(int j=i+1; j<dim; ++j)
        {
            EXPECT_DOUBLE_EQ(0.0, X(i,j));
            EXPECT_DOUBLE_EQ(0.0, X(j,i));
        }
    }

}


TEST(Utilities_LyapunovSolver, DampedOscillator)
{
    const int dim = 2;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim,dim);
    A << 0.0, 1.0,
        -1.0, -0.1;
    
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(dim,dim);
    

    LyapunovSolver<double> solver;
    solver.compute(A, Q);

    auto& X = solver.matrixX();

    
    Eigen::MatrixXcd resid = A.cast<std::complex<double>>().adjoint()*X + X*A.cast<std::complex<double>>() + Q.cast<std::complex<double>>();

    for(int i=0; i<dim; ++i){
        for(int j=0; j<dim; ++j){
            EXPECT_NEAR(0.0, resid(i,j).real(),6e-15);
        }
    }
}

TEST(Utilities_LyapunovSolver, MaternTest)
{
    const int dim = 2;
    
    Eigen::MatrixXd A(2,2);
    A << 0,            1.0,
        -133.333333, -23.094;

    Eigen::MatrixXd Q(2,2);
    Q << 0.0, 0.0,
         0.0, 6158.4;

    LyapunovSolver<double> solver;
    solver.compute(A, Q);

    auto& X = solver.matrixX();

    
    Eigen::MatrixXcd resid = A.cast<std::complex<double>>().adjoint()*X + X*A.cast<std::complex<double>>() + Q.cast<std::complex<double>>();

    for(int i=0; i<dim; ++i){
        for(int j=0; j<dim; ++j){
            EXPECT_NEAR(0.0, resid(i,j).real(),1e-10);
        }
    }
}

