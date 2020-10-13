/**
OVERVIEW:
The goal of this example is to fit a monotone function to stress-strain
data obtained during the 2017 Sandia Fracture Challenge.  For our area of
interest, the stress should be a monotonically increasing function of the
strain.

The MonotoneExpansion class provides a way of characterizing
monotone functions and can be fit to data in a least squares sense using
the Gauss-Newton algorithm, which is implemented in the "FitData" function
below.
*/

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"

#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/HDF5/H5Object.h"

#include <Eigen/Dense>

using namespace muq::Approximation;
using namespace muq::Utilities;

/** Reads SFC3 stress strain data. */
std::pair<Eigen::VectorXd, Eigen::VectorXd> ReadData()
{
  auto f = muq::Utilities::OpenFile("data/LTA01.h5");

  std::cout << "Reading strain data..." << std::endl;
  std::cout << "  Units: " << std::string(f["/Strain"].attrs["Units"]) << std::endl;
  std::cout << "  Size:  " << f["/Strain"].rows() << std::endl;

  Eigen::VectorXd strain = f["/Strain"];

  std::cout << "Reading stress data..." << std::endl;
  std::cout << "  Units: " << std::string(f["/Stress"].attrs["Units"]) << std::endl;
  std::cout << "  Size:  " << f["/Stress"].rows() << std::endl;

  Eigen::VectorXd stress = f["/Stress"];

  return std::make_pair(strain, stress);
};


std::shared_ptr<MonotoneExpansion> SetupExpansion(unsigned order)
{
    auto poly = std::make_shared<Legendre>();
    std::vector<std::shared_ptr<IndexedScalarBasis>> bases(1, poly);
    std::shared_ptr<MultiIndexSet> multis = MultiIndexFactory::CreateTotalOrder(1, order);

    Eigen::MatrixXd polyCoeffs = Eigen::MatrixXd::Zero(1,multis->Size());
    polyCoeffs(0,1) = 500.0; // add an initial linear trend

    auto polyBase = std::make_shared<BasisExpansion>(bases, multis, polyCoeffs);
    auto expansion = std::make_shared<MonotoneExpansion>(polyBase);

    return expansion;
}

/** Find parameters of the montone expansion that minimize the L2 norm between
    the predictions and observations, i.e., solve the nonlinear least squares
    problem for the parameters describing the monotone function.
*/
void FitData(Eigen::VectorXd             const& x,
             Eigen::VectorXd             const& y,
             std::shared_ptr<MonotoneExpansion> expansion)
{
  Eigen::VectorXd coeffs = expansion->GetCoeffs();

  Eigen::VectorXd preds(x.size());
  Eigen::VectorXd newPreds(x.size());
  Eigen::VectorXd newCoeffs;

  Eigen::VectorXd resid, newResid, step, xslice;
  Eigen::MatrixXd jac(x.size(), coeffs.size());

  const int maxLineIts = 10;
  const int maxIts = 20;

  // Use the current monotone parameterization to make predictions at every point
  for(int k=0; k<x.size(); ++k){
    xslice = x.segment(k,1);
    preds(k) = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(xslice,coeffs).at(0))(0);
  }

  resid = y-preds;
  double sse = resid.squaredNorm();

  for(int i=0; i<maxIts; ++i){

    // Compute the jacobian at the current point
    for(int k=0; k<x.size(); ++k){
      xslice = x.segment(k,1);
      jac.row(k) = boost::any_cast<Eigen::MatrixXd>(expansion->Jacobian(1,0,xslice,coeffs));
    }

    // Compute the Gauss-Newton step
    step = jac.colPivHouseholderQr().solve(resid);
    newCoeffs = coeffs + step;

    // Use the current monotone parameterization to make predictions at every point
    for(int k=0; k<x.size(); ++k){
      xslice = x.segment(k,1);
      newPreds(k) = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(xslice,newCoeffs).at(0))(0);
    }

    newResid = y-newPreds;
    double newsse = newResid.squaredNorm();

    // Backtracing line search to guarantee a sufficient descent
    int lineIt = 0;
    while((newsse > sse - 1e-7)&&(lineIt<maxLineIts)){
      step *= 0.5;
      newCoeffs = coeffs + step;

      // Compute the residuals at the new point
      for(int k=0; k<x.size(); ++k){
        xslice = x.segment(k,1);
        newPreds(k) = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(xslice,newCoeffs).at(0))(0);
      }

      newResid = y-newPreds;
      newsse = newResid.squaredNorm();
    }

    if(lineIt == maxLineIts){
      std::cout << "WARNING: Line search failed, terminating Gauss-Newton optimizer." << std::endl;
      return;
    }

    // The line search was successful, so update the coefficients and residuals
    coeffs = newCoeffs;
    preds = newPreds;
    sse = newsse;
    resid = newResid;

    std::cout << "Iteration " << i << ", SSE = "<< sse << std::endl;
  }
}

void WriteResults(Eigen::VectorXd             const& x,
                  std::shared_ptr<MonotoneExpansion> expansion)
{
  Eigen::VectorXd preds(x.size());
  Eigen::VectorXd xslice;

  for(int k=0; k<x.size(); ++k){
    xslice = x.segment(k,1);
    preds(k) = boost::any_cast<Eigen::VectorXd>(expansion->Evaluate(xslice).at(0))(0);
  }

  auto f = muq::Utilities::OpenFile("results/StressPredictions.h5");
  f["/Strain"] = x;
  f["/Stress"] = preds;
}


int main()
{

  Eigen::VectorXd strain, stress;
  std::tie(strain, stress) = ReadData();

  unsigned polyOrder = 7;
  auto expansion = SetupExpansion(polyOrder);

  FitData(strain, stress, expansion);

  WriteResults(strain, expansion);

  return 0;
};
