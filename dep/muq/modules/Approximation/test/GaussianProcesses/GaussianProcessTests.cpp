#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <iostream>

using namespace muq::Approximation;
using namespace muq::Utilities;
using namespace muq::Modeling;

TEST(Approximation_GP, Discretize)
{
  auto kernel = SquaredExpKernel(1, 2.0, 0.35);

  // Create the GP
  ZeroMean mean(1, 1);
  GaussianProcess gp(mean, kernel);

  Eigen::RowVectorXd pts = Eigen::RowVectorXd::LinSpaced(100, 0, 1);

  std::shared_ptr<Gaussian> gauss = gp.Discretize(pts);

  Eigen::MatrixXd gaussCov = gauss->GetCovariance();
  Eigen::VectorXd gaussMean = gauss->GetMean();

  Eigen::MatrixXd gpMean, gpCov;
  std::tie(gpMean, gpCov) = gp.Predict(pts, GaussianProcess::FullCov);

  for(int i=0; i<gpMean.size(); ++i){
    EXPECT_DOUBLE_EQ(gpMean(i), gaussMean(i));
    for(int j=0; j<=i; ++j)
      EXPECT_DOUBLE_EQ(gpCov(i,j), gaussCov(i,j));
  }

}


TEST(Approximation_GP, DiscretizeConcatenate)
{
  auto kernel1 = SquaredExpKernel(1, 2.0, 0.35);
  auto kernel2 = MaternKernel(1, 2.0, 0.2, 5.0/2.0);
  auto kernel = Concatenate(kernel1, kernel2);

  // Create the GP
  ZeroMean mean(2, 2);
  GaussianProcess gp(mean, kernel);

  Eigen::RowVectorXd pts = Eigen::RowVectorXd::LinSpaced(100, 0, 1);

  std::shared_ptr<Gaussian> gauss = gp.Discretize(pts);

  Eigen::MatrixXd gaussCov = gauss->GetCovariance();
  Eigen::VectorXd gaussMean = gauss->GetMean();

  Eigen::MatrixXd gpMean, gpCov;
  std::tie(gpMean, gpCov) = gp.Predict(pts, GaussianProcess::FullCov);

  for(int i=0; i<gpMean.size(); ++i){
    EXPECT_DOUBLE_EQ(gpMean(i), gaussMean(i));
    for(int j=0; j<=i; ++j)
      EXPECT_DOUBLE_EQ(gpCov(i,j), gaussCov(i,j));
  }

}


TEST(Approximation_GP, HyperFit1d)
{
    const unsigned numPred  = 50;
    const unsigned maxTrain = 10;

    const double pi = 4.0 * atan(1.0);

    // Set linearly space locations where we want to evaluate the Gaussian Process
    Eigen::RowVectorXd predLocs(numPred);
    predLocs.setLinSpaced(numPred, 0,1);

    // Generate random training locations
    Eigen::RowVectorXd trainLocs = Eigen::RowVectorXd::LinSpaced(maxTrain, 0, 1);

    // Generate the training data (with some random noise)
    Eigen::RowVectorXd trainData(maxTrain);
    for(int i=0; i<maxTrain; ++i)
	    trainData(i) = sin(4*2*pi*trainLocs(i) );

    trainData += sqrt(1e-4)*RandomGenerator::GetNormal(maxTrain).transpose();

    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35, {0.1,10}, {0, 100} ) * PeriodicKernel(dim, 1.0, 0.75, 0.25, {0.5,5.0}, {0.5,5.0}, {0.25,0.5}) + WhiteNoiseKernel(dim, 1e-3, {1e-8,100});

    // Create the GP
    ZeroMean mean(dim, 1);
    GaussianProcess gp(mean, kernel);


    for(int i=0; i<trainData.size(); ++i)
        gp.Condition(trainLocs.col(i), trainData.col(i));

    // Fit the hyperparameters
    //gp.Optimize();

    // Make a prediction
    Eigen::MatrixXd postMean = gp.PredictMean(predLocs);

    //for(int i=0; i<predLocs.cols(); ++i)
    //	std::cout << "[" << predLocs(0,i) << ", " << postMean(0,i) << "]," << std::endl;

}

TEST(Approximation_GP, HyperFit2d)
{
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    const unsigned numPred  = 50;
    const unsigned maxTrain = 50;

    // Set linearly space locations where we want to evaluate the Gaussian Process
    Eigen::MatrixXd predLocs(2, numPred);
    for(int i=0; i<numPred; ++i)
    {
	    predLocs(0,i) = uniform_dist(e1);
	    predLocs(1,i) = uniform_dist(e1);
    }

    // Generate random training locations
    Eigen::MatrixXd trainLocs(2, maxTrain);

    for(int i=0; i<maxTrain; ++i)
    {
	    trainLocs(0,i) = uniform_dist(e1);
	    trainLocs(1,i) = uniform_dist(e1);
    }

    // Generate the training data (with some random noise)
    Eigen::RowVectorXd trainData(maxTrain);
    for(int i=0; i<maxTrain; ++i)
	    trainData(i) = trainLocs.col(i).squaredNorm();

    // define a tensor product kernel
    std::vector<unsigned> inds1{0};
    std::vector<unsigned> inds2{1};
    const unsigned dim = 2;
    auto kernel = SquaredExpKernel(dim, inds1, 2.0, 0.35, {0.1,10} , {0.0, 100})*SquaredExpKernel(dim, inds2, 2.0, 0.35, {0.1,10}, {0.0, 100.0} );

    // Create the GP
    ZeroMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);

    for(int i=0; i<trainLocs.cols(); ++i)
        gp.Condition(trainLocs.col(i),trainData.col(i));

    //gp.Optimize();

    // Make a prediction
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> post = gp.Predict(predLocs, GaussianProcess::FullCov);

}

TEST(Approximation_GP, ValueCondition)
{

    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35);

    const double delta =1e-4;
    Eigen::MatrixXd evalLocs(dim, 2);
    evalLocs << 0, delta;

    // Create the GP
    ZeroMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);

    Eigen::VectorXd loc(1);
    loc << 0.0;

    Eigen::VectorXd val(1);
    val << 1.0;

    Eigen::MatrixXd valCov(1,1);
    valCov(0,0) = 0.0;

    auto H = std::make_shared<IdentityOperator>(dim);

    auto obs = std::make_shared<ObservationInformation>(H, loc, val, valCov);

    gp.Condition(obs);

    Eigen::MatrixXd field = gp.PredictMean(evalLocs);

    EXPECT_NEAR(val(0), field(0,0), 1e-6);
}


TEST(Approximation_GP, GradientCondition)
{
    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35);

    const double delta = 1e-4;
    Eigen::MatrixXd evalLocs(dim, 2);
    evalLocs << 0, delta;

    // Create the GP
    ZeroMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);

    Eigen::VectorXd loc(1);
    loc << 0.0;

    Eigen::VectorXd derivVal(1);
    derivVal << 0.5;

    Eigen::MatrixXd derivCov(1,1);
    derivCov(0,0) = 1e-4;

    std::vector<std::vector<int>> derivCoords;
    derivCoords.push_back({0});

    auto H = std::make_shared<IdentityOperator>(dim);

    auto derivObs = std::make_shared<DerivativeObservation>(H, loc, derivVal, derivCov, derivCoords);

    Eigen::MatrixXd crossCov(1,1);
    derivObs->FillCrossCov(evalLocs.col(1),kernel.Clone(),crossCov);
    EXPECT_DOUBLE_EQ(kernel.GetPosDerivative(evalLocs.col(1), evalLocs.col(0), {0})(0,0), crossCov(0,0));

    gp.Condition(derivObs);

    Eigen::MatrixXd field = gp.PredictMean(evalLocs);
    EXPECT_NEAR( derivVal(0), (field(1)-field(0))/delta, 1e-5);

    // field = gp.Sample(evalLocs);
    // EXPECT_NEAR( derivVal(0), (field(1)-field(0))/delta, 1e-5);

    Eigen::MatrixXd cov;
    std::tie(field, cov) = gp.Predict(evalLocs,GaussianProcess::FullCov);
    EXPECT_NEAR( derivVal(0), (field(1)-field(0))/(evalLocs(0,1)-evalLocs(0,0)), 1e-5);
}


TEST(Approximation_GP, HessianCondition)
{
    const unsigned dim = 1;
    auto kernel = SquaredExpKernel(dim, 2.0, 0.35);

    const double delta = 1e-4;
    Eigen::MatrixXd evalLocs(dim, 4);
    evalLocs << 0, delta, 2.0*delta, 3.0*delta;

    // Create the GP
    ZeroMean mean(dim, 1);
    auto gp = ConstructGP(mean, kernel);

    Eigen::VectorXd loc(1);
    loc << 0.0;

    // Value of the second derivative
    Eigen::VectorXd derivVal(1);
    derivVal << 0.5;

    // Observation variance
    Eigen::MatrixXd derivCov(1,1);
    derivCov(0,0) = 1e-4;

    std::vector<std::vector<int>> derivCoords;
    derivCoords.push_back({0,0});

    auto H = std::make_shared<IdentityOperator>(dim);

    auto derivObs = std::make_shared<DerivativeObservation>(H, loc, derivVal, derivCov, derivCoords);

    Eigen::MatrixXd crossCov(1,1);
    derivObs->FillCrossCov(evalLocs.col(1),kernel.Clone(),crossCov);
    EXPECT_DOUBLE_EQ(kernel.GetPosDerivative(evalLocs.col(1), evalLocs.col(0), {0,0})(0,0), crossCov(0,0));

    gp.Condition(derivObs);

    Eigen::MatrixXd field = gp.PredictMean(evalLocs);
    double fdDeriv = (2.0*field(0) - 5.0*field(1) + 4.0*field(2) -field(3))/(delta*delta);
    EXPECT_NEAR( derivVal(0), fdDeriv, 1e-5);

//    field = gp.Sample(evalLocs);
//    fdDeriv = (2.0*field(0) - 5.0*field(1) + 4.0*field(2) -field(3))/(delta*delta);
//    EXPECT_NEAR( derivVal(0), fdDeriv, 1e-5);

    Eigen::MatrixXd cov, field2;
    std::tie(field2, cov) = gp.Predict(evalLocs,GaussianProcess::FullCov);
    fdDeriv = (2.0*field2(0) - 5.0*field2(1) + 4.0*field2(2) -field2(3))/(delta*delta);
    EXPECT_NEAR( derivVal(0), fdDeriv, 1e-5);

}
