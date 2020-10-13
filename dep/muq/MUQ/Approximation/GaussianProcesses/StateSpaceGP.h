#ifndef STATESPACEGP_H
#define STATESPACEGP_H

#include <memory>

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearSDE.h"

#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"

#include <Eigen/Core>

namespace muq
{
namespace Approximation
{


class StateSpaceGP : public GaussianProcess
{
public:

    StateSpaceGP(MeanFunctionBase&           meanIn,
                 KernelBase&                 kernelIn,
                 boost::property_tree::ptree options = boost::property_tree::ptree()) : StateSpaceGP(meanIn.Clone(), kernelIn.Clone(), options){};

    StateSpaceGP(std::shared_ptr<MeanFunctionBase> meanIn,
                 std::shared_ptr<KernelBase>       covKernelIn,
                 boost::property_tree::ptree       options = boost::property_tree::ptree());

    virtual ~StateSpaceGP() = default;
    
    virtual Eigen::MatrixXd Sample(Eigen::MatrixXd const& times) override;


    virtual std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Predict(Eigen::MatrixXd const& newLocs,
                                                                CovarianceType         covType) override;

    virtual Eigen::MatrixXd PredictMean(Eigen::MatrixXd const& newPts) override;


    virtual double LogLikelihood(Eigen::MatrixXd const& xs,
                                 Eigen::MatrixXd const& vals) override;

    virtual double MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                         bool                        computeGrad = true) override;


    std::shared_ptr<muq::Modeling::LinearSDE> GetSDE(){return sde;};

    std::shared_ptr<muq::Modeling::LinearOperator> GetObs(){return obsOp;};

    void SetObs(std::shared_ptr<muq::Modeling::LinearOperator> newObs);

    Eigen::MatrixXd GetCov(){return L.triangularView<Eigen::Lower>()*L.transpose();};



    const int stateDim;

    //static std::shared_ptr<StateSpaceGP> Concatenate(std::vector<std::shared_ptr<StateSpaceGP>> const& gps);

private:

    StateSpaceGP(std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> ssInfo,
                 std::shared_ptr<MeanFunctionBase> meanIn,
                 std::shared_ptr<KernelBase>       covKernelIn);


    void SortObservations();

    bool ComputeAQ(double dt);
    Eigen::MatrixXd sdeA, sdeQ;
    double dtAQ; // the last deltat passed to the ComputeAQ function

    // Stocastic Differential equation describing correlations
    std::shared_ptr<muq::Modeling::LinearSDE> sde;

    // Observation operator
    std::shared_ptr<muq::Modeling::LinearOperator> obsOp;

    // Cholesky factor of the stationary covariance matrix (used to initialize SDE integration)
    Eigen::MatrixXd L;

    std::shared_ptr<MeanFunctionBase> mean;
    std::shared_ptr<KernelBase>       covKernel;


};


} // namespace Approximation
} // namespace GP





#endif
