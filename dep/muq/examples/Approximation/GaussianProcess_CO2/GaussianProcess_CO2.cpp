#include <iostream>

#include "MUQ/Utilities/HDF5/H5Object.h"

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"
#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"
#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

using namespace std;
using namespace muq::Utilities;
using namespace muq::Approximation;


template<typename MeanType, typename KernelType>
std::shared_ptr<GaussianProcess> BuildGP(MeanType const& mu, KernelType const& kernel)
{
    return std::make_shared<GaussianProcess>(mu.Clone(), kernel.Clone());
}

template<typename MeanType, typename KernelType>
std::shared_ptr<GaussianProcess> BuildStateSpaceGP(MeanType const& mu, KernelType const& kernel)
{
    boost::property_tree::ptree options;
    options.put("SDE.dt", 1e-2);

    return std::make_shared<StateSpaceGP>(mu.Clone(), kernel.Clone(), options);
}

int main()
{
    // Open and read in data
    string dataFile = "data/MaunaLoaCO2.h5";
    H5Object f = OpenFile(dataFile);

    Eigen::VectorXd times          = f["/Weekly/Dates" ];
    Eigen::VectorXd concentrations = f["/Weekly/Concentrations" ];

    // Long term trend
    double k1_var = 360;
    double k1_length = 60;
    double k1_nu = 5.0/2.0;
    auto k1 = MaternKernel(1, k1_var, k1_length, k1_nu);

    // Periodic component
    double k2_var = 1.0;
    double k2_period = 1.0;
    double k2_length = 0.5;
    auto k2 = PeriodicKernel(1, k2_var, k2_length, k2_period);

    // Quasi-periodic term
    double k3_var = 10.0;
    double k3_length = 30;
    double k3_nu = 3.0/2.0;
    auto k3 = MaternKernel(1, k3_var, k3_length, k3_nu);

    // Short term trends or anomalis
    double k4_var = 1.0;
    double k4_length = 2.0;
    double k4_nu = 3.0/2.0;
    auto k4 = MaternKernel(1, k4_var, k4_length, k4_nu);

    // Combine the individual kernels
    auto k = k1 + k2*k3 + k4;

    // Add a linear mean function
    LinearMean mu(1.8, 370-2000.0*1.8);

    // Construct a Gaussian Process using standard Matrix backend
    //auto gp = BuildGP(mu, k);

    // Construct a Gaussian Processing using a StateSpace formulation
    auto gp = BuildStateSpaceGP(mu,k);

    // Define pediction points
    int numPts = 800;
    Eigen::MatrixXd evalPts(1, numPts);
    evalPts.row(0) = Eigen::VectorXd::LinSpaced(numPts, 2010, 2020);

    gp->Condition(times.transpose(), concentrations.transpose(), 16);
    Eigen::MatrixXd postMean, postVar;
    std::tie(postMean,postVar) = gp->Predict(evalPts, GaussianProcess::DiagonalCov);

    // Write results to new file
    string writeFile = "results/CO2_Prediction.h5";
    H5Object fout = OpenFile(writeFile);

    fout["/Predict/Dates"] = evalPts;
    fout["/Predict/Concentrations"] = postMean;
    fout["/Predict/ConcentrationVariance"] = postVar;
    
    return 0;
}
