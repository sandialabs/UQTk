#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;


TEST(Approximation_GP, LinearTransformKernel)
{

    const unsigned dim = 2;
    Eigen::MatrixXd sigma2(2,2);
    sigma2 << 1.0, 0.9,
	            0.9, 1.5;

    auto kernel = ConstantKernel(dim, sigma2) * SquaredExpKernel(dim, 2.0, 0.35 );

    EXPECT_EQ(2, kernel.coDim);

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3,2);

    auto kernel2 = A*kernel;
    auto kernel3 = A * ConstantKernel(dim, sigma2) * SquaredExpKernel(dim, 2.0, 0.35 );

    Eigen::VectorXd x1(dim);
    x1 << 0.1, 0.4;

    Eigen::VectorXd x2(dim);
    x2 << 0.2, 0.7;

    Eigen::MatrixXd result2 = kernel2.Evaluate(x1,x2);
    Eigen::MatrixXd result3 = kernel3.Evaluate(x1,x2);

    Eigen::MatrixXd expected = A * kernel.Evaluate(x1,x2) * A.transpose();

    for(int j=0; j<A.rows(); ++j)
    {
	    for(int i=0; i<A.rows(); ++i)
	    {
	      EXPECT_NEAR(expected(i,j), result2(i,j), 1e-15);
	      EXPECT_NEAR(expected(i,j), result3(i,j), 1e-15);
	    }
    }
}



TEST(Approximation_GP, Clone)
{

    const unsigned dim = 2;
    auto kernel = ConstantKernel(dim, 2.0) * SquaredExpKernel(dim, 2.0, 0.35 );

    std::shared_ptr<KernelBase> kernel_copy = kernel.Clone();


    EXPECT_DOUBLE_EQ(kernel.inputDim, kernel_copy->inputDim);
    EXPECT_DOUBLE_EQ(kernel.coDim, kernel_copy->coDim);
    EXPECT_DOUBLE_EQ(kernel.numParams, kernel_copy->numParams);


    Eigen::VectorXd x1(dim);
    x1 << 0.1, 0.4;

    Eigen::VectorXd x2(dim);
    x2 << 0.2, 0.7;

    Eigen::MatrixXd result = kernel.Evaluate(x1,x2);
    Eigen::MatrixXd result_ptr = kernel_copy->Evaluate(x1,x2);

    EXPECT_DOUBLE_EQ(result(0,0), result_ptr(0,0));
}

TEST(Approximation_GP, KernelConcatenation)
{
    const double sigma2 = 2.0;
    const double length = 0.5;

    MaternKernel kernel12(1, sigma2, length, 1.0/2.0);
    MaternKernel kernel32(1, sigma2, length, 3.0/2.0);

    auto kernel2 = Concatenate(kernel12, kernel32);

}


// TEST(Approximation_GP, SeperableProduct)
// {
//     {

//         const unsigned dim = 2;
//         auto kernel1 = SquaredExpKernel(dim, 2.0, 3.5) * SquaredExpKernel(dim, 1.0, 0.5);

//         auto comps1 = kernel1.GetSeperableComponents();
//         EXPECT_EQ(1, comps1.size());
//     }

//     {
//         std::vector<unsigned> inds1{0};
//         std::vector<unsigned> inds2{1};
//         const unsigned dim = 2;
//         auto kernel2 = SquaredExpKernel(dim, inds1, 2.0, 0.35, {0.1,10} ) * SquaredExpKernel(dim, inds2, 2.0, 0.35, {0.1,10} );

//         auto comps2 = kernel2.GetSeperableComponents();
//         EXPECT_EQ(2, comps2.size());
//     }

//     {
//         std::vector<unsigned> inds1{0};
//         std::vector<unsigned> inds2{1,2};

//         const unsigned dim = 3;
//         auto kernel2 = SquaredExpKernel(dim, inds1, 2.0, 0.35, {0.1,10} ) * SquaredExpKernel(dim, inds2, 2.0, 0.35, {0.1,10} );

//         auto comps2 = kernel2.GetSeperableComponents();
//         EXPECT_EQ(2, comps2.size());
//     }

//     {
//         std::vector<unsigned> inds1{0,1};
//         std::vector<unsigned> inds2{1,2};

//         const unsigned dim = 3;
//         auto kernel2 = SquaredExpKernel(dim, inds1, 2.0, 0.35, {0.1,10} ) * SquaredExpKernel(dim, inds2, 2.0, 0.35, {0.1,10} );

//         auto comps2 = kernel2.GetSeperableComponents();
//         EXPECT_EQ(1, comps2.size());
//     }

// }


// TEST(Approximation_GP, StateSpaceError)
// {
//     std::vector<std::shared_ptr<KernelBase>> kernels;
//     kernels.push_back( std::make_shared<SquaredExpKernel>(1, 1.0, 1.0) );
//     kernels.push_back( std::make_shared<SquaredExpKernel>(2, 1.0, 1.0) );
//
//     CoregionalKernel kernel(2, Eigen::MatrixXd::Identity(2,2), kernels);
//
//     EXPECT_THROW(kernel.GetStateSpace(), muq::NotImplementedError);
//
// }

TEST(Approximation_GP, SquaredExpKernel)
{
    Eigen::VectorXd pt1(1);
    pt1 << 0.5;

    Eigen::VectorXd pt2(1);
    pt2 << 0.75;

    const double sigma2 = 2.0;
    const double length = 0.5;

    auto kernel = SquaredExpKernel(1, sigma2, length);

    Eigen::MatrixXd cov = kernel.Evaluate(pt1,pt1);
    EXPECT_DOUBLE_EQ(sigma2, cov(0,0));

    cov = kernel.Evaluate(pt1,pt2);
    EXPECT_DOUBLE_EQ(sigma2*exp(-0.5 * std::pow((pt1(0)-pt2(0))/length,2.0)), cov(0,0));


    // Finite difference derivative test
    const double eps = 1e-5;
    Eigen::VectorXd leftPt = pt1;
    leftPt(0) -= eps;
    Eigen::VectorXd rightPt = pt1;
    rightPt(0) += eps;

    double centerVal = kernel.Evaluate(pt1,pt2)(0,0);
    double leftVal = kernel.Evaluate(leftPt,pt2)(0,0);
    double rightVal = kernel.Evaluate(rightPt,pt2)(0,0);

    // Check first derivative
    double deriv = kernel.GetPosDerivative(pt1,pt2,{0})(0,0);
    EXPECT_NEAR((rightVal-leftVal)/(2.0*eps), deriv, 1e-4);

    // Check second derivative
    deriv = kernel.GetPosDerivative(pt1,pt2,{0,0})(0,0);
    EXPECT_NEAR((rightVal - 2.0*centerVal + leftVal)/(eps*eps), deriv, 1e-4);
}


double MaternKernel12(double r, double l)
{
  return exp(-r/l);
}

double MaternKernel32(double r, double l)
{
  return (1.0 + sqrt(3.0)*r/l)*exp(-sqrt(3)*r/l);
}

double MaternKernel52(double r, double l)
{
  return (1.0 + sqrt(5.0)*r/l + 5.0*r*r/(3.0*l*l))*exp(-sqrt(5)*r/l);
}


TEST(Approximation_GP, MaternKernel)
{
    Eigen::VectorXd pt1(1);
    pt1 << 0.5;

    Eigen::VectorXd pt2(1);
    pt2 << 0.5;

    const double sigma2 = 2.0;
    const double length = 0.5;

    EXPECT_THROW(MaternKernel(1, sigma2, length, 2.0), std::invalid_argument);

    {
      MaternKernel kernel12(1, sigma2, length, 1.0/2.0);
      double truth = MaternKernel12((pt2-pt1).norm(), length);
      EXPECT_NEAR(sigma2*truth, kernel12.Evaluate(pt1,pt2)(0,0), 1e-10);

      MaternKernel kernel32(1, sigma2, length, 3.0/2.0);
      truth = MaternKernel32((pt2-pt1).norm(), length);
      EXPECT_NEAR(sigma2*truth, kernel32.Evaluate(pt1,pt2)(0,0), 1e-10);

      MaternKernel kernel52(1, sigma2, length, 5.0/2.0);
      truth = MaternKernel52((pt2-pt1).norm(), length);
      EXPECT_NEAR(sigma2*truth, kernel52.Evaluate(pt1,pt2)(0,0), 1e-10);
    }

    pt2 << 0.75;

    MaternKernel kernel12(1, sigma2, length, 1.0/2.0);
    double truth = MaternKernel12((pt2-pt1).norm(), length);
    EXPECT_NEAR(sigma2*truth, kernel12.Evaluate(pt1,pt2)(0,0), 1e-10);

    MaternKernel kernel32(1, sigma2, length, 3.0/2.0);
    truth = MaternKernel32((pt2-pt1).norm(), length);
    EXPECT_NEAR(sigma2*truth, kernel32.Evaluate(pt1,pt2)(0,0), 1e-10);

    MaternKernel kernel52(1, sigma2, length, 5.0/2.0);
    truth = MaternKernel52((pt2-pt1).norm(), length);
    EXPECT_NEAR(sigma2*truth, kernel52.Evaluate(pt1,pt2)(0,0), 1e-10);


    // Finite difference derivative test
    const double eps = 1e-6;
    Eigen::VectorXd leftPt = pt1;
    leftPt(0) -= eps;
    Eigen::VectorXd rightPt = pt1;
    rightPt(0) += eps;

    double centerVal = kernel52.Evaluate(pt1,pt2)(0,0);
    double leftVal = kernel52.Evaluate(leftPt,pt2)(0,0);
    double rightVal = kernel52.Evaluate(rightPt,pt2)(0,0);

    // Check first derivative
    double deriv = kernel52.GetPosDerivative(pt1,pt2,{0})(0,0);
    EXPECT_NEAR((rightVal-leftVal)/(2.0*eps), deriv, 1e-4);

    // Check second derivative
    deriv = kernel52.GetPosDerivative(pt1,pt2,{0,0})(0,0);
    EXPECT_NEAR((rightVal - 2.0*centerVal + leftVal)/(eps*eps), deriv, 1e-4);

}
