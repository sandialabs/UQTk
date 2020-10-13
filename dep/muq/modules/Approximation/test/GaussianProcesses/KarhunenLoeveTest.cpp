
#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Utilities;
using namespace muq::Approximation;

TEST(Approximation_GP, KarhunenLoeve_BasicConstruction)
{

    double sigma2 = 1.0;
    double length = 0.2;
    double nu = 3.0/2.0;

    auto kernel = std::make_shared<MaternKernel>(1, sigma2, length, nu);

    const int numSeeds = 300;
    const int numEval = 200;

    const double lb = 0.0;
    const double ub = 1.0;

    Eigen::MatrixXd seedPts(1,numSeeds);
    for(int i=0; i<numSeeds; ++i)
        seedPts(0, i) = lb + (ub-lb)*double(i)/double(numSeeds-1);

    Eigen::VectorXd seedWts = (ub-lb)/double(numSeeds)*Eigen::VectorXd::Ones(numSeeds);
    seedWts(0) *= 0.5;
    seedWts(numSeeds-1) *= 0.5;

    // Construct the Karhunen-Loeve decomposition based on a discrete eigenvalue problem at the seed points
    boost::property_tree::ptree options;
    options.put("KarhunenLoeve.TruncationType", "FixedNumber");
    options.put("KarhunenLoeve.EnergyTol", 0.999);
    options.put("KarhunenLoeve.NumModes", 100);
    KarhunenLoeveExpansion kl(kernel, seedPts, seedWts, options);

    // Evaluate the KL expansion at a bunch of other points
    Eigen::MatrixXd evalPts(1,numEval);
    for(int i=0; i<numEval; ++i)
        evalPts(0, i) = lb + (ub-lb)*double(i)/double(numEval-1);

    Eigen::MatrixXd modes = kl.GetModes(evalPts);

    const unsigned int N = 1e4;

    //double var = (samps-Eigen::VectorXd::Constant(N,mu)).array().pow(2.0).sum()/(N-1.0);
    Eigen::VectorXd vars = modes.array().square().rowwise().sum();
    for(int i=0; i<numEval; ++i)
      EXPECT_NEAR(vars(i), sigma2, 1.0e-3);

    // check the sample mean
    Eigen::VectorXd samps(N);
    const double x = 0.25;
    for( unsigned int i=0; i<N; ++i )
      samps(i) = kl.Evaluate(Eigen::VectorXd::Constant(1, x), RandomGenerator::GetNormal(modes.cols()))(0);

    double mu = samps.mean();
    EXPECT_NEAR(mu, 0.0, 3.0e-2);

    // check the variance
    double var = (samps-Eigen::VectorXd::Constant(N,mu)).array().pow(2.0).sum()/(N-1.0);

    EXPECT_NEAR(var, sigma2, 3.0e-2);
}


TEST(Approximation_GP, KarhunenLoeve_GaussQuad)
{
    const double lb = 0.0;
    const double ub = 1.0;

    // Get Gauss-Legendre points on [lb,ub]
    GaussQuadrature gq(std::make_shared<Legendre>());
    gq.Compute(50);
    Eigen::VectorXd gaussPts = gq.Points().transpose();
    gaussPts = 0.5*(gaussPts+Eigen::VectorXd::Ones(gaussPts.size()));
    gaussPts = (ub-lb)*gaussPts + lb*Eigen::VectorXd::Ones(gaussPts.size());

    Eigen::VectorXd gaussWts = 0.5*(ub-lb)*gq.Weights();

    double sigma2 = 1.0;
    double length = 0.2;
    double nu = 3.0/2.0;

    auto kernel = std::make_shared<MaternKernel>(1, sigma2, length, nu);

    // Construct the Karhunen-Loeve decomposition based on a discrete eigenvalue problem at the seed points
    KarhunenLoeveExpansion kl(kernel, gaussPts.transpose(), gaussWts);

    // Evaluate the KL expansion at a bunch of other points
    const int numEval = 200;
    Eigen::MatrixXd evalPts(1,numEval);
    for(int i=0; i<numEval; ++i)
        evalPts(0, i) = lb + (ub-lb)*double(i)/double(numEval-1);

    Eigen::MatrixXd modes = kl.GetModes(evalPts);

    Eigen::MatrixXd klCov = modes*modes.transpose();
    Eigen::MatrixXd trueCov = kernel->BuildCovariance(evalPts);

    for(int j=0; j<numEval; ++j){
      for(int i=0; i<numEval; ++i){
        EXPECT_NEAR(trueCov(i,j), klCov(i,j), 1e-2);
      }
    }
}
