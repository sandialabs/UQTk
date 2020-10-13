#include <gtest/gtest.h>

#include <Eigen/Core>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/ModPiece.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

TEST(GaussianDistributionTests, SpecifyMean) {

  // create the distributions
  auto standard1D = std::make_shared<Gaussian>(1);

  const Eigen::VectorXd mu = Eigen::VectorXd::Ones(3);
  auto standard = std::make_shared<Gaussian>(mu);
  EXPECT_EQ(standard1D->Dimension(), 1);
  EXPECT_EQ(standard->Dimension(), 3);

  // compute the log density for a standard normal
  Eigen::VectorXd x(1);
  x << 2.5;
  Eigen::VectorXd x3(3);
  x3 << 2.0, 1.4, 3.0;

  EXPECT_DOUBLE_EQ(standard1D->LogDensity(x), -0.5*std::log(2.0*M_PI) - 0.5*x(0)*x(0));
  EXPECT_DOUBLE_EQ(standard->LogDensity(x3), -0.5*mu.size()*std::log(2.0*M_PI)-0.5*(x3-mu).dot(x3-mu));

  const unsigned int N = 1.0e5;

  Eigen::VectorXd mean1D = standard1D->Sample();
  Eigen::VectorXd mean = standard->Sample();
  for( unsigned int i=0; i<N; ++i ) {
    mean1D += standard1D->Sample();
    mean += standard->Sample();
  }
  mean1D /= (N+1.0);
  mean /= (N+1.0);

  EXPECT_NEAR(mean1D(0), 0.0, 1.0e-2);
  EXPECT_NEAR((mean-mu).norm(), 0.0, 1.0e-2);
}

TEST(GaussianDistributionTests, ChangeHyperparameters) {

  // create a standard normal (one dimension)
  auto standard = std::make_shared<Gaussian>(1, Gaussian::Mean);

  const unsigned int N = 1.0e5;

  // compute the log density for a standard normal
  Eigen::VectorXd x(1);
  x << 2.5;

  Eigen::VectorXd mu(1);
  mu << 1.0;

  Eigen::MatrixXd cov(1,1);
  cov << 0.2;

  EXPECT_DOUBLE_EQ(standard->LogDensity(x, mu), -0.5*std::log(2.0*M_PI)-(mu(0)-x(0))*(mu(0)-x(0))/2.0);

  Eigen::VectorXd samps(N);
  for( unsigned int i=0; i<N; ++i ) {
    samps(i) = standard->Sample(mu)(0);
  }

  double mean = samps.mean();
  double sigma2 = ((samps.array()-mean)*(samps.array()-mean)).sum()/(N-1);
  EXPECT_NEAR(mean, 1.0, 8.0e-2);
  EXPECT_NEAR(sigma2, 1.0, 8.0e-2);

  standard = std::make_shared<Gaussian>(1, Gaussian::Mean | Gaussian::FullCovariance);
  Eigen::VectorXd covVec = Eigen::Map<Eigen::VectorXd>(cov.data(), cov.rows()*cov.cols()).eval();
  EXPECT_DOUBLE_EQ(standard->LogDensity(x, mu, covVec), -0.5*(std::log(2.0*M_PI)+std::log(cov(0,0)))-(mu(0)-x(0))*(mu(0)-x(0))/(2.0*cov(0,0)));

  for( unsigned int i=0; i<N; ++i ) {
    samps(i) = standard->Sample(mu, covVec)(0);
  }
  mean = samps.sum()/(N+1.0);
  samps = samps.array()-mean;
  sigma2 = (samps.array()*samps.array()).sum()/N;
  EXPECT_NEAR(mean, 1.0, 8.0e-2);
  EXPECT_NEAR(sigma2, 0.25, 8.0e-2);

  // EXPECT_DOUBLE_EQ(standard->LogDensity(x, prec), -0.5*(std::log(2.0*M_PI)+std::log(0.5))-(1.0-x(0))*(1.0-x(0)));
  // samps(0) = boost::any_cast<double>(standard->Sample(prec));
  // for( unsigned int i=0; i<N; ++i ) {
  //   samps(i+1) = boost::any_cast<double>(standard->Sample(prec));
  // }
  // mean = samps.sum()/(N+1.0);
  // samps = samps.array()-mean;
  // sigma2 = (samps.array()*samps.array()).sum()/N;
  // EXPECT_NEAR(mean, 1.0, 1.0e-2);
  // EXPECT_NEAR(sigma2, 0.5, 1.0e-2);
}

TEST(GaussianDistributionTests, SpecifyBoth) {
  // create the distributions
  const Eigen::Vector3d cov = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector3d prec = Eigen::Vector3d::Random().cwiseAbs();
  const Eigen::Vector3d mu = Eigen::Vector3d::Ones();
  auto covDist = std::make_shared<Gaussian>(mu, cov);
  auto precDist = std::make_shared<Gaussian>(mu, prec, Gaussian::Mode::Precision);
  EXPECT_EQ(covDist->Dimension(), 3);
  EXPECT_EQ(precDist->Dimension(), 3);

  // compute the log density for a standard normal
  Eigen::VectorXd x(3);
  x << 2.0, 1.4, 3.0;
  EXPECT_DOUBLE_EQ(covDist->LogDensity(x), -0.5*(3.0*std::log(2.0*M_PI)+cov.array().log().sum())-(x-mu).dot((1.0/cov.array()).matrix().asDiagonal()*(x-mu))/2.0);
  EXPECT_DOUBLE_EQ(precDist->LogDensity(x), -0.5*(3.0*std::log(2.0*M_PI)-prec.array().log().sum())-(x-mu).dot(prec.asDiagonal()*(x-mu))/2.0);

  const unsigned int N = 8.0e5;

  Eigen::VectorXd meanCov = boost::any_cast<Eigen::VectorXd>(covDist->Sample());
  Eigen::VectorXd meanPrec = boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  for( unsigned int i=0; i<N; ++i ) {
    meanCov += boost::any_cast<Eigen::VectorXd>(covDist->Sample());
    meanPrec += boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  }
  meanCov /= (N+1.0);
  meanPrec /= (N+1.0);

  EXPECT_NEAR((meanCov-mu).norm(), 0.0, 1.0e-2);
  EXPECT_NEAR((meanPrec-mu).norm(), 0.0, 1.0e-2);
}


TEST(GaussianDistributionTests, ModPieceInputs) {

  int dim = 2;
  auto dens = std::make_shared<Gaussian>(dim, Gaussian::Mean | Gaussian::FullCovariance)->AsDensity();

  EXPECT_EQ(3, dens->inputSizes.size());
  EXPECT_EQ(dim, dens->inputSizes(0));
  EXPECT_EQ(dim, dens->inputSizes(1));
  EXPECT_EQ(dim*dim, dens->inputSizes(2));

  Eigen::VectorXd xIn = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd muIn = Eigen::VectorXd::Ones(dim);
  Eigen::VectorXd covIn = Eigen::VectorXd::Ones(dim*dim);

  double logDens = dens->Evaluate(xIn, muIn, covIn).at(0)(0);

}

TEST(GaussianDistributionTests, Gradient) {

  int dim = 20;

  Eigen::VectorXd mu = Eigen::VectorXd::Random(dim);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(dim,dim);
  Eigen::MatrixXd cov = A*A.transpose() + 1e-8*Eigen::MatrixXd::Identity(dim,dim);

  std::shared_ptr<Distribution> dist = std::make_shared<Gaussian>(mu,cov);

  Eigen::VectorXd testPt = Eigen::VectorXd::Ones(dim);

  Eigen::VectorXd grad = dist->GradLogDensity(0, testPt);
  Eigen::VectorXd trueGrad = -1.0*cov.llt().solve(testPt-mu);

  for(int i=0; i<dim; ++i)
    EXPECT_NEAR(trueGrad(i), grad(i),5e-12);
}


TEST(GaussianDistributionTests, DiagonalCovPrec) {
  // the covariance or precision scale (variance in the cov. case and 1/variance in the prec. case)
  Eigen::VectorXd covDiag = 1e-2 * Eigen::VectorXd::Ones(2) + RandomGenerator::GetUniform(2);
  Eigen::VectorXd precDiag = covDiag.array().inverse();

  Eigen::VectorXd mu = Eigen::VectorXd::Zero(2);
  auto diagCov = std::make_shared<Gaussian>(mu, covDiag, Gaussian::Mode::Covariance);
  auto diagPrec = std::make_shared<Gaussian>(mu, precDiag, Gaussian::Mode::Precision);
  EXPECT_EQ(diagCov->Dimension(), 2);
  EXPECT_EQ(diagPrec->Dimension(), 2);

  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  EXPECT_DOUBLE_EQ(diagCov->LogDensity(x), -0.5*(2.0*std::log(2.0*M_PI)+covDiag.array().log().sum())-x.dot((1.0/covDiag.array()).matrix().asDiagonal()*x)/2.0);
  EXPECT_DOUBLE_EQ(diagPrec->LogDensity(x), -0.5*(2.0*std::log(2.0*M_PI)-precDiag.array().log().sum())-x.dot(precDiag.asDiagonal()*x)/2.0);

  const unsigned int N = 1.0e5;

  RandomGenerator::SetSeed(5192012);
  Eigen::VectorXd meanCov = diagCov->Sample();

  RandomGenerator::SetSeed(5192012);
  Eigen::VectorXd meanPrec = diagCov->Sample();

  EXPECT_NEAR(meanCov(0),meanPrec(0),1e-12);
  EXPECT_NEAR(meanCov(1),meanPrec(1),1e-12);

  for( unsigned int i=0; i<N-1; ++i ) {
    meanCov += diagCov->Sample();
    meanPrec += diagPrec->Sample();
  }
  meanCov /= N;
  meanPrec /= N;

  EXPECT_NEAR(0.0, meanCov(0), 1e-2);
  EXPECT_NEAR(0.0, meanCov(1), 1e-2);
  EXPECT_NEAR(0.0, meanPrec(0), 1e-2);
  EXPECT_NEAR(0.0, meanPrec(1), 1e-2);
}

TEST(GaussianDistributionTests, MatrixCovPrec) {

  const unsigned int dim = 5;
  Eigen::MatrixXd cov = Eigen::MatrixXd::Random(dim, dim);
  cov = Eigen::MatrixXd::Identity(dim, dim) + cov*cov.transpose();

  Eigen::MatrixXd prec = Eigen::MatrixXd::Random(dim, dim);
  prec = Eigen::MatrixXd::Identity(dim, dim) + prec*prec.transpose();

  Eigen::LLT<Eigen::MatrixXd> covChol(cov);
  const Eigen::MatrixXd covL = covChol.matrixL();

  Eigen::LLT<Eigen::MatrixXd> precChol(prec);
  const Eigen::MatrixXd precL = precChol.matrixL();

  Eigen::VectorXd mu = Eigen::VectorXd::Zero(dim);
  auto covDist = std::make_shared<Gaussian>(mu, cov, Gaussian::Mode::Covariance);
  auto precDist = std::make_shared<Gaussian>(mu, prec, Gaussian::Mode::Precision);
  EXPECT_EQ(covDist->Dimension(), dim);
  EXPECT_EQ(Gaussian::Mode::Covariance, covDist->GetMode());
  EXPECT_EQ(Gaussian::Mode::Precision, precDist->GetMode());

  const Eigen::VectorXd x = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd delta = covL.triangularView<Eigen::Lower>().solve(x);
  covL.triangularView<Eigen::Lower>().transpose().solveInPlace(delta);
  EXPECT_DOUBLE_EQ(covDist->LogDensity(x), -0.5*(dim*std::log(2.0*M_PI)+2.0*covL.diagonal().array().log().sum())-x.dot(delta)/2.0);
  EXPECT_DOUBLE_EQ(precDist->LogDensity(x), -0.5*(dim*std::log(2.0*M_PI)-2.0*precL.diagonal().array().log().sum())-x.dot(prec*x)/2.0);


  const unsigned int N = 5.0e5;
  Eigen::VectorXd meanCov = boost::any_cast<Eigen::VectorXd>(covDist->Sample());
  Eigen::VectorXd meanPrec = boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  for( unsigned int i=0; i<N-1; ++i ) {
    meanCov += boost::any_cast<Eigen::VectorXd>(covDist->Sample());
    meanPrec += boost::any_cast<Eigen::VectorXd>(precDist->Sample());
  }
  meanCov /= N;
  meanPrec /= N;

  EXPECT_NEAR(0.0, meanCov(0), 1e-2);
  EXPECT_NEAR(0.0, meanCov(1), 1e-2);
  EXPECT_NEAR(0.0, meanPrec(0), 1e-2);
  EXPECT_NEAR(0.0, meanPrec(1), 1e-2);
}
