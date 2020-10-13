#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"
#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/ProductOperator.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Inference/Filtering/KalmanFilter.h"
#include "MUQ/Inference/Filtering/KalmanSmoother.h"

#include "MUQ/Utilities/Exceptions.h"

#include <Eigen/Dense>

using namespace muq::Approximation;
using namespace muq::Modeling;
using namespace muq::Utilities;
using namespace muq::Inference;

StateSpaceGP::StateSpaceGP(std::shared_ptr<MeanFunctionBase> meanIn,
                           std::shared_ptr<KernelBase>       covKernelIn,
                           boost::property_tree::ptree       options) : StateSpaceGP(covKernelIn->GetStateSpace(options),
                                                                                     meanIn,
                                                                                     covKernelIn)
{
};


StateSpaceGP::StateSpaceGP(std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> ssInfo,
                           std::shared_ptr<MeanFunctionBase> meanIn,
                           std::shared_ptr<KernelBase>       covKernelIn) : GaussianProcess(meanIn, covKernelIn), stateDim(std::get<2>(ssInfo).rows()),
                                                                            sde(std::get<0>(ssInfo)),
                                                                            obsOp(std::get<1>(ssInfo)),
                                                                            L(std::get<2>(ssInfo).selfadjointView<Eigen::Lower>().llt().matrixL()),
                                                                            mean(meanIn),
                                                                            covKernel(covKernelIn)
{
}

void StateSpaceGP::SortObservations()
{
    if(hasNewObs)
    {
        auto obsComp = [](std::shared_ptr<ObservationInformation> a, std::shared_ptr<ObservationInformation> b) { return a->loc(0) < b->loc(0); };
        std::sort(observations.begin(), observations.end(), obsComp);
        hasNewObs = false;
    }
}

Eigen::MatrixXd StateSpaceGP::Sample(Eigen::MatrixXd const& times)
{

    // Make space for the simulated GP
    Eigen::MatrixXd output(obsOp->rows(), times.size());

    if(observations.size()==0){

        // Generate sample for initial condition
        Eigen::VectorXd x = L.triangularView<Eigen::Lower>()*RandomGenerator::GetNormal(L.rows());

        output.col(0) = obsOp->Apply(x);

        // Step through the each time and integrate the SDE between times
        for(int i=0; i<times.size()-1; ++i)
        {
            x = sde->EvolveState(x, times(i+1)-times(i));
            output.col(i+1) = obsOp->Apply(x);
        }

    }else{

        throw muq::NotImplementedError("The Sample function of muq::Approximation::StateSpaceGP does not currently support Gaussian Processes that have been conditioned on data.");
    }

    return output;
}


Eigen::MatrixXd StateSpaceGP::PredictMean(Eigen::MatrixXd const& newPts)
{
    return Predict(newPts, GaussianProcess::NoCov).first;
}

bool StateSpaceGP::ComputeAQ(double dt)
{
    bool computeAQ = false;
    if( (sdeA.rows()==0) || (sdeA.cols()==0) || (std::abs(dt-dtAQ)>2.0*std::numeric_limits<double>::epsilon()) ){
        computeAQ = true;
    }

    if(computeAQ){
        std::tie(sdeA, sdeQ) = sde->Discretize(dt);
        dtAQ = dt;
    }

    return computeAQ;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> StateSpaceGP::Predict(Eigen::MatrixXd const& times,
                                                                  CovarianceType         covType)
{

    if(observations.size()==0)
        return GaussianProcess::Predict(times, covType);


    // Make sure the evaluations points are one dimensional
    if(times.rows() != 1)
        throw muq::WrongSizeError("In StateSpaceGP::Predict: The StateSpaceGP class only supports 1d fields, but newPts has " + std::to_string(times.rows()) + " dimensions (i.e., rows).");

    // Make sure the evaluation points are sorted
    for(int i=1; i<times.cols(); ++i)
    {
        if(times(0,i)<=times(0,i-1))
            throw std::invalid_argument("In StateSpaceGP::Predict: The input points must be monotonically increasing, but index " + std::to_string(i-1) + " has a value of " + std::to_string(times(0,i-1)) + " and index " + std::to_string(i) + " has a value of " + std::to_string(times(1,i)) + ".");
    }

    if(covType==GaussianProcess::FullCov)
    {
        throw std::invalid_argument("In StateSpaceGP::Predict: The statespace GP does not compute the full covariance.");
    }

    std::vector< std::pair<Eigen::VectorXd, Eigen::MatrixXd>> evalDists(times.cols());

    std::vector< std::pair<Eigen::VectorXd, Eigen::MatrixXd>> obsDists(observations.size());
    std::vector< std::pair<Eigen::VectorXd, Eigen::MatrixXd>> obsFilterDists(observations.size());

    std::pair<Eigen::VectorXd, Eigen::MatrixXd>* currDist;
    double currTime, nextTime;

    int obsInd = 0;
    int evalInd = 0;


    // Is the first time an observation or an evaluation?
    if(observations.at(0)->loc(0) < times(0)){
        currTime = observations.at(0)->loc(0);

        obsDists.at(0).second = L.triangularView<Eigen::Lower>()*L.transpose();
        obsDists.at(0).first = Eigen::VectorXd::Zero(stateDim);

        auto H = std::make_shared<ProductOperator>(observations.at(0)->H, obsOp);
        obsFilterDists.at(0) = KalmanFilter::Analyze(obsDists.at(0),
                                                     H,
                                                     observations.at(0)->obs - mean->Evaluate(observations.at(0)->loc),
                                                     observations.at(0)->obsCov);

        currDist = &obsFilterDists.at(0);

        obsInd++;
    }else{
        currTime = times(0);

        evalDists.at(0).second = L.triangularView<Eigen::Lower>()*L.transpose();
        evalDists.at(0).first = Eigen::VectorXd::Zero(stateDim);

        currDist = &evalDists.at(0);

        evalInd++;
    }

    // Loop through the remaining times
    for(int i=1; i<times.cols() + observations.size(); ++i)
    {
        bool hasObs = false;
        if(obsInd<observations.size()){
            if(evalInd>=times.cols()){
                hasObs = true;
            }else if(observations.at(obsInd)->loc(0) < times(evalInd)){
                hasObs = true;
            }
        }

        // Is this time an observation or an evaluation?
        if(hasObs){

            ComputeAQ(observations.at(obsInd)->loc(0) - currTime);

            obsDists.at(obsInd).first = sdeA * currDist->first;
            obsDists.at(obsInd).second = sdeA * currDist->second * sdeA.transpose() + sdeQ; //= sde->EvolveDistribution(*currDist, observations.at(obsInd)->loc(0) - currTime);

            currTime = observations.at(obsInd)->loc(0);

            auto H = std::make_shared<ProductOperator>(observations.at(obsInd)->H, obsOp);
            obsFilterDists.at(obsInd) = KalmanFilter::Analyze(obsDists.at(obsInd),
                                                              H,
                                                              observations.at(obsInd)->obs - mean->Evaluate(observations.at(obsInd)->loc),
                                                              observations.at(obsInd)->obsCov);

            currDist = &obsFilterDists.at(obsInd);

            obsInd++;

        }else{

            ComputeAQ(times(evalInd) - currTime);

            evalDists.at(evalInd).first = sdeA * currDist->first;
            evalDists.at(evalInd).second = sdeA * currDist->second * sdeA.transpose() + sdeQ;

            currTime = times(evalInd);
            currDist = &evalDists.at(evalInd);

            evalInd++;
        }
    }

    ///////////////////////////////////////////////////
    // BACKWARD (SMOOTHING) PASS
    ///////////////////////////////////////////////////

    // Find the last evaluation point before (or at the same time) as the last observation
    obsInd = observations.size()-2;
    evalInd = times.size()-1;
    for(evalInd = times.size()-1; evalInd >=0; --evalInd)
    {
        if(times(evalInd) <= observations.at(observations.size()-1)->loc(0))
            break;
    }

    if(evalInd>=0)
    {
        nextTime = observations.at(observations.size()-1)->loc(0);

        std::pair<Eigen::VectorXd, Eigen::MatrixXd> filterDist = obsDists.at(observations.size()-1);
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> smoothDist = obsFilterDists.at(observations.size()-1);

        while( evalInd>=0 )
        {
            // Have we passed any new observations on the way backward?
            bool hasObs = false;
            if(obsInd>=0){
                if(observations.at(obsInd)->loc(0) > times(evalInd) )
                    hasObs = true;
            }

            if( hasObs )
            {
                // Update the value at the observation point if it's not the last observation

                ComputeAQ(nextTime - observations.at(obsInd)->loc(0));
                smoothDist = KalmanSmoother::Analyze(obsFilterDists.at(obsInd), filterDist, smoothDist, sdeA);
                filterDist.swap( obsDists.at(obsInd) );

                nextTime = observations.at(obsInd)->loc(0);
                obsInd--;

            }else{
                ComputeAQ(nextTime - times(evalInd));

                smoothDist = KalmanSmoother::Analyze(evalDists.at(evalInd), filterDist, smoothDist, sdeA);

                filterDist.swap( evalDists.at(evalInd) );
                evalDists.at(evalInd) = smoothDist;

                nextTime = times(evalInd);
                evalInd--;
            }
        }
    }

    //////////////////////////////////////////////////
    // GET MEAN AND COVARIANCE FROM SDE SOLUTION
    //////////////////////////////////////////////////

    // Copy the solution into the output
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> output;
    output.first.resize(coDim, times.cols());
    for(int i=0; i<times.cols(); ++i)
        output.first.col(i) = mean->Evaluate(times.col(i)) + obsOp->Apply(evalDists.at(i).first);

    if(covType==GaussianProcess::BlockCov){

        output.second.resize(coDim, coDim*times.cols());

        for(int i=0; i<times.cols(); ++i)
            output.second.block(0,i*coDim,coDim,coDim) = obsOp->Apply( obsOp->Apply(evalDists.at(i).second).transpose() );

    }else if(covType==GaussianProcess::DiagonalCov){
        output.second.resize(coDim, times.cols());

        for(int i=0; i<times.cols(); ++i)
            output.second.col(i) = obsOp->Apply( obsOp->Apply(evalDists.at(i).second).transpose() ).diagonal();
    }

    return output;
}


double StateSpaceGP::LogLikelihood(Eigen::MatrixXd const& xs,
                                   Eigen::MatrixXd const& vals)
{
    return 0.0;
}

double StateSpaceGP::MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                           bool                        computeGrad)
{
    return 0.0;
}





void StateSpaceGP::SetObs(std::shared_ptr<muq::Modeling::LinearOperator> newObs)
{

    if(newObs->cols() != obsOp->cols())
        throw muq::WrongSizeError("In StateSpaceGP::SetObs: The new observation operator has " + std::to_string(newObs->cols()) + " columns, which does not match the system dimension " + std::to_string(obsOp->cols()));

    obsOp = newObs;
}

// std::shared_ptr<StateSpaceGP> StateSpaceGP::Concatenate(std::vector<std::shared_ptr<StateSpaceGP>> const& gps)
// {
//     // Build the concatenated SDE
//     std::vector<std::shared_ptr<muq::Modeling::LinearSDE>> sdes(gps.size());
//     for(int i=0; i<gps.size(); ++i)
//         sdes.at(i) = gps.at(i)->GetSDE();

//     auto sde = LinearSDE::Concatenate(sdes);

//     // Build a concatenated observation operator
//     std::vector<std::shared_ptr<muq::Modeling::LinearOperator>> obsOps(gps.size());
//     for(int i=0; i<gps.size(); ++i)
//         obsOps.at(i) = gps.at(i)->GetObs();

//     auto H = std::make_shared<BlockDiagonalOperator>(obsOps);

//     // Build the concatenated covariance
//     Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(sde->stateDim, sde->stateDim);
//     int currRow = 0;
//     for(int i=0; i<gps.size(); ++i)
//     {
//         Q.block(currRow,currRow, gps.at(i)->stateDim, gps.at(i)->stateDim) = gps.at(i)->GetCov();
//         currRow += gps.at(i)->stateDim;
//     }

//     return std::make_shared<StateSpaceGP>(sde, H, Q);
// }
