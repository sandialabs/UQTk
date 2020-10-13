#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include "MUQ/Utilities/Exceptions.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>


using namespace muq::Approximation;


TEST(Approximation_GP, MaternStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.15;

    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    ZeroMean mu(1,1);
    StateSpaceGP gp(mu, kernel);

    EXPECT_EQ(nu+0.5, gp.stateDim);

    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(100, 0, 1);
    Eigen::MatrixXd realization = gp.Sample(obsTimes);

}

TEST(Approximation_GP, StateSpace_DistributionIntegration)
{

    const double sigma2 = 1.0;
    const double length = 0.15;

    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    std::shared_ptr<muq::Modeling::LinearSDE> sde; // The underying SDE
    std::shared_ptr<muq::Modeling::LinearOperator> H; // Observation matrix to go from SDE state to actual GP variable
    Eigen::MatrixXd pinf; // Steady state covariance matrix

    std::tie(sde, H, pinf) = kernel.GetStateSpace();

    Eigen::MatrixXd p0 = Eigen::MatrixXd::Identity(sde->stateDim, sde->stateDim);
    Eigen::MatrixXd mu0 = Eigen::VectorXd::Ones(sde->stateDim);

    // Integrate the SDE for a while
    Eigen::MatrixXd muT, pT;
    std::tie(muT, pT) = sde->EvolveDistribution(mu0,p0,30);

    EXPECT_NEAR(0.0, muT(0), 1e-14);
    EXPECT_NEAR(0.0, muT(1), 1e-14);

    EXPECT_NEAR(pinf(0,0), pT(0,0), 1e-11);
    EXPECT_NEAR(pinf(0,1), pT(0,1), 1e-11);
    EXPECT_NEAR(pinf(1,0), pT(1,0), 1e-11);
    EXPECT_NEAR(pinf(1,1), pT(1,1), 1e-11);
}

TEST(Approximation_GP, StateSpace_DistributionIntegration2)
{

    const double sigma2 = 1.0;
    const double length = 0.15;

    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    std::shared_ptr<muq::Modeling::LinearSDE> sde; // The underying SDE
    std::shared_ptr<muq::Modeling::LinearOperator> H; // Observation matrix to go from SDE state to actual GP variable
    Eigen::MatrixXd pinf; // Steady state covariance matrix

    std::tie(sde, H, pinf) = kernel.GetStateSpace();

    Eigen::MatrixXd p0 = Eigen::MatrixXd::Identity(sde->stateDim, sde->stateDim);
    Eigen::VectorXd mu0 = Eigen::VectorXd::Ones(sde->stateDim);
    std::pair<Eigen::VectorXd, Eigen::MatrixXd> dist0 = std::make_pair(mu0,p0);

    // Integrate the SDE for a while
    Eigen::MatrixXd muT, pT;
    std::tie(muT, pT) = sde->EvolveDistribution(dist0,30);

    EXPECT_NEAR(0.0, muT(0), 1e-14);
    EXPECT_NEAR(0.0, muT(1), 1e-14);

    EXPECT_NEAR(pinf(0,0), pT(0,0), 1e-11);
    EXPECT_NEAR(pinf(0,1), pT(0,1), 1e-11);
    EXPECT_NEAR(pinf(1,0), pT(1,0), 1e-11);
    EXPECT_NEAR(pinf(1,1), pT(1,1), 1e-11);
}

TEST(Approximation_GP, PeriodicStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.6;
    const double period = 0.25;
    const double periodN = 50; // how many steps per period

    PeriodicKernel kernel(1, sigma2, length, period);

    boost::property_tree::ptree options;
    options.put("PeriodicKernel.StateSpace.NumTerms",8);
    options.put("SDE.dt", 5e-5);

    ZeroMean mu(1,1);
    StateSpaceGP gp(mu, kernel, options);

    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(5*periodN+1, 0, 5*period);

    Eigen::MatrixXd realization = gp.Sample(obsTimes);

    // Make sure the sample is periodic
    for(int i=0; i<obsTimes.size()-periodN-1; ++i)
        EXPECT_NEAR(realization(0,i), realization(0,i+periodN), 1e-1);

}

TEST(Approximation_GP, SumStateSpace)
{

    MaternKernel kernel1(1, 3.0, 1.0, 7.0/2.0);
    MaternKernel kernel2(1, 0.15, 0.25, 1.0/2.0);

    auto kernel = kernel1 + kernel2;

    boost::property_tree::ptree options;
    options.put("PeriodicKernel.StateSpace.NumTerms",8);
    options.put("SDE.dt", 5e-5);

    ZeroMean mu(1,1);
    StateSpaceGP gp(mu, kernel, options);

    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(100, 0, 4);

    Eigen::MatrixXd realization = gp.Sample(obsTimes);
}

TEST(Approximation_GP, ProductStateSpace)
{

    const double sigma2 = 1.0;
    const double length = 0.6;
    const double nu = 3.0/2.0;
    const double period = 0.25;
    const double periodN = 50; // how many steps per period

    PeriodicKernel kernel1(1, sigma2, 0.8, period);
    MaternKernel kernel2(1, sigma2, 2.0, nu);

    auto kernel12 = kernel1*kernel2;
    auto kernel21 = kernel2*kernel1;
    auto kernel22 = kernel2*kernel2;

    boost::property_tree::ptree options;
    options.put("PeriodicKernel.StateSpace.NumTerms",7);
    options.put("SDE.dt", 1e-4);

    ZeroMean mu(1,1);
    StateSpaceGP gp1(mu, kernel1, options);
    StateSpaceGP gp2(mu, kernel2, options);
    EXPECT_EQ(int(nu+0.5), gp2.stateDim);

    StateSpaceGP gp12(mu, kernel12, options);
    StateSpaceGP gp21(mu, kernel21, options);

    EXPECT_THROW(kernel22.GetStateSpace(options), muq::NotImplementedError);
    EXPECT_EQ(gp1.stateDim *gp2.stateDim, gp12.stateDim);


    // draw a random sample from the SDE model
    Eigen::VectorXd obsTimes = Eigen::VectorXd::LinSpaced(5*periodN+1, 0, 5*period);

    Eigen::MatrixXd realization1  = gp1.Sample(obsTimes);
    Eigen::MatrixXd realization2  = gp2.Sample(obsTimes);
    Eigen::MatrixXd realization12 = gp12.Sample(obsTimes);

}

TEST(Approximation_GP, StateSpacePredict_Interpolation)
{

    const double sigma2 = 1.0;
    const double length = 0.3;
    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-5);

    ZeroMean mu(1,1);
    StateSpaceGP gp1(mu, kernel, options);
    GaussianProcess gp2(mu, kernel);

    const int numEvals = 100;
    Eigen::MatrixXd evalPts(1,numEvals);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numEvals, 0, 2);

    // Make a prediction about the prior
    Eigen::MatrixXd predMu, predCov;
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);

    Eigen::MatrixXd predMu2, predCov2;
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    // Condition both GPs with some data
    Eigen::MatrixXd obsLoc = 0.5*Eigen::MatrixXd::Ones(1,1);
    Eigen::MatrixXd obsData = 0.25*Eigen::MatrixXd::Ones(1,1);
    const double obsVar = 1e-2;

    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // add another observation
    obsLoc(0) = 1.5;
    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // Make a posterior prediction
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    const double meanTol = 1e-2;
    const double covTol = 1e-2;
    for(int j=0; j<predMu.size(); ++j){
        EXPECT_NEAR(predMu2(0,j), predMu(0,j), meanTol);
        EXPECT_NEAR(predCov2(0,j), predCov(0,j), covTol);
    }
}

TEST(Approximation_GP, StateSpacePredict_ExtrapolateRight)
{

    const double sigma2 = 1.0;
    const double length = 0.3;
    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-5);

    ZeroMean mu(1,1);
    StateSpaceGP gp1(mu, kernel, options);
    GaussianProcess gp2(mu, kernel);

    const int numEvals = 100;
    Eigen::MatrixXd evalPts(1,numEvals);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numEvals, 1.6, 3);

    // Make a prediction about the prior
    Eigen::MatrixXd predMu, predCov;
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);

    Eigen::MatrixXd predMu2, predCov2;
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    // Condition both GPs with some data
    Eigen::MatrixXd obsLoc = 0.5*Eigen::MatrixXd::Ones(1,1);
    Eigen::MatrixXd obsData = 0.25*Eigen::MatrixXd::Ones(1,1);
    const double obsVar = 1e-2;

    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // add another observation
    obsLoc(0) = 1.5;
    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // Make a posterior prediction
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    const double meanTol = 1e-2;
    const double covTol = 1e-2;
    for(int j=0; j<predMu.size(); ++j){
        EXPECT_NEAR(predMu2(0,j), predMu(0,j), meanTol);
        EXPECT_NEAR(predCov2(0,j), predCov(0,j), covTol);
    }

}

TEST(Approximation_GP, StateSpacePredict_ExtrapolateLeft)
{

    const double sigma2 = 1.0;
    const double length = 0.3;
    const double nu = 3.0/2.0;

    MaternKernel kernel(1, sigma2, length, nu);

    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-5);

    ZeroMean mu(1,1);
    StateSpaceGP gp1(mu, kernel, options);
    GaussianProcess gp2(mu, kernel);

    const int numEvals = 100;
    Eigen::MatrixXd evalPts(1,numEvals);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numEvals, -1, 0.2);

    // Make a prediction about the prior
    Eigen::MatrixXd predMu, predCov;
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);

    Eigen::MatrixXd predMu2, predCov2;
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    // Condition both GPs with some data
    Eigen::MatrixXd obsLoc = 0.5*Eigen::MatrixXd::Ones(1,1);
    Eigen::MatrixXd obsData = 0.25*Eigen::MatrixXd::Ones(1,1);
    const double obsVar = 1e-2;

    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // add another observation
    obsLoc(0) = 1.5;
    gp1.Condition(obsLoc, obsData, obsVar);
    gp2.Condition(obsLoc, obsData, obsVar);

    // Make a posterior prediction
    std::tie(predMu, predCov) = gp1.Predict(evalPts, GaussianProcess::DiagonalCov);
    std::tie(predMu2, predCov2) = gp2.Predict(evalPts, GaussianProcess::DiagonalCov);

    const double meanTol = 1e-2;
    const double covTol = 1e-2;
    for(int j=0; j<predMu.size(); ++j){
        EXPECT_NEAR(predMu2(0,j), predMu(0,j), meanTol);
        EXPECT_NEAR(predCov2(0,j), predCov(0,j), covTol);
    }
}
