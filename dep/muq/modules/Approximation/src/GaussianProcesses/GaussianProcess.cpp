#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Approximation;
using namespace muq::Modeling;
using namespace muq::Utilities;

double muq::Approximation::nlopt_obj(unsigned n, const double *x, double *nlopt_grad, void *opt_info)
{

    Eigen::Map<const Eigen::VectorXd> params(x,n);
    double logLikely;

    OptInfo* info = (OptInfo*) opt_info;

    info->gp->Kernel()->SetParams(params);

    if(nlopt_grad){
        Eigen::Map<Eigen::VectorXd> grad(nlopt_grad, n);
        logLikely = info->gp->MarginalLogLikelihood(grad);
    }else{
        logLikely = info->gp->MarginalLogLikelihood();
    }

    return logLikely;
}


GaussianProcess::GaussianProcess(std::shared_ptr<MeanFunctionBase> meanIn,
				 std::shared_ptr<KernelBase>       covKernelIn) :
    mean(meanIn),
    covKernel(covKernelIn),
    inputDim(covKernel->inputDim),
    coDim(covKernel->coDim),
    hasNewObs(false)
{};

GaussianProcess& GaussianProcess::Condition(Eigen::Ref<const Eigen::MatrixXd> const& loc,
                                            Eigen::Ref<const Eigen::MatrixXd> const& vals,
                                            double                                   obsVar)
{
    auto H = std::make_shared<IdentityOperator>(coDim);

    assert(loc.cols()>0);
    assert(vals.cols()==loc.cols());
    assert(vals.rows()==coDim);
    assert(obsVar>=0);

    for(int i=0; i<loc.cols(); ++i){
        auto obs = std::make_shared<ObservationInformation>(H,
                                                            loc.col(i),
                                                            vals.col(i),
                                                            obsVar*Eigen::MatrixXd::Identity(coDim,coDim));

        Condition(obs);
    }

    return *this;
}

GaussianProcess& GaussianProcess::Condition(std::shared_ptr<ObservationInformation> obs)
{
    observations.push_back(obs);
    hasNewObs = true;

    return *this;
}

std::shared_ptr<muq::Modeling::Gaussian> GaussianProcess::Discretize(Eigen::MatrixXd const& pts)
{
  Eigen::MatrixXd mean;
  Eigen::MatrixXd cov;

  std::tie(mean,cov) =  Predict(pts, GaussianProcess::FullCov);

  Eigen::Map<Eigen::VectorXd> meanMap(mean.data(), mean.rows()*mean.cols());

  return std::make_shared<muq::Modeling::Gaussian>(meanMap,cov);
}


void GaussianProcess::ProcessObservations()
{
    if((hasNewObs)&&(observations.size()>0)){

        obsDim = 0;
        for(int i=0; i<observations.size(); ++i)
            obsDim += observations.at(i)->H->rows();

        // Build the covariance between observation locations
        Eigen::MatrixXd trainCov(obsDim, obsDim);
        int currCol = 0;
        for(int j=0; j<observations.size(); ++j)
        {
            int currRow = 0;
            for(int i=0; i<j; ++i){
                observations.at(i)->FillCrossCov(observations.at(j), covKernel, trainCov.block(currRow, currCol, observations.at(i)->H->rows(), observations.at(j)->H->rows()));
                currRow += observations.at(i)->H->rows();
            }
            //covKernel->BuildCovariance(observations.at(j)->loc, observations.at(j)->loc);
            observations.at(j)->FillSelfCov(covKernel, trainCov.block(currRow, currRow, observations.at(j)->H->rows(), observations.at(j)->H->rows()));
            currCol += observations.at(j)->H->rows();
        }

        trainCov.triangularView<Eigen::Lower>() = trainCov.triangularView<Eigen::Upper>().transpose();

        covSolver = trainCov.selfadjointView<Eigen::Lower>().ldlt();

        // Evaluate the mean function
        Eigen::VectorXd trainDiff(obsDim);
        int currRow = 0;
        for(int i=0; i<observations.size(); ++i){
            auto derivObs = std::dynamic_pointer_cast<DerivativeObservation>(observations.at(i));
            if(derivObs){
              trainDiff.segment(currRow, observations.at(i)->H->rows()) = observations.at(i)->obs - observations.at(i)->H->Apply(mean->GetDerivative(observations.at(i)->loc, derivObs->derivCoords));
            }else{
              trainDiff.segment(currRow, observations.at(i)->H->rows()) = observations.at(i)->obs - observations.at(i)->H->Apply(mean->Evaluate(observations.at(i)->loc));
            }
            currRow += observations.at(i)->H->rows();
        }

        sigmaTrainDiff = covSolver.solve(trainDiff);

        hasNewObs = false;
    }
}


void GaussianProcess::Optimize()
{
    ProcessObservations();

    OptInfo info;
    info.gp = this;

    Eigen::VectorXd params = covKernel->GetParams();

    const Eigen::MatrixXd bounds = covKernel->GetParamBounds();
    Eigen::VectorXd lbs = bounds.row(0);
    Eigen::VectorXd ubs = bounds.row(1);


    // nlopt_opt opt;

    // opt = nlopt_create(NLOPT_LD_LBFGS, params.rows()); /* algorithm and dimensionality */
    // nlopt_set_lower_bounds(opt, lbs.data());
    // nlopt_set_upper_bounds(opt, ubs.data());
    // nlopt_set_vector_storage(opt, 100);
    // nlopt_set_max_objective(opt, muq::Approximation::nlopt_obj, (void*) &info);

    // double maxLikely;
    // nlopt_optimize(opt, params.data(), &maxLikely);

    // covKernel->SetParams(params);

    // nlopt_destroy(opt);
}


Eigen::MatrixXd GaussianProcess::BuildCrossCov(Eigen::MatrixXd const& newLocs)
{
    Eigen::MatrixXd crossCov( obsDim, coDim*newLocs.cols());
    int currCol = 0;

    for(int j=0; j<observations.size(); ++j){
        for(int i=0; i<newLocs.cols(); ++i){
          observations.at(j)->FillCrossCov(newLocs.col(i), covKernel, crossCov.block(currCol, i*coDim, observations.at(j)->H->rows(), coDim));

            // Eigen::MatrixXd temp = covKernel
            // Eigen::MatrixXd temp = covKernel->BuildCovariance(newLocs.col(i), observations.at(j)->loc);
            // Eigen::MatrixXd temp2 = observations.at(j)->H->Apply( covKernel->BuildCovariance(newLocs.col(i), observations.at(j)->loc));
            // crossCov.block(i*coDim, currCol, coDim, observations.at(j)->H->rows()) = observations.at(j)->H->Apply( covKernel->BuildCovariance(newLocs.col(i), observations.at(j)->loc)).transpose();
        }
        currCol += observations.at(j)->H->rows();
    }

    return crossCov.transpose();
}
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> GaussianProcess::Predict(Eigen::MatrixXd const& newLocs,
                                                                     CovarianceType         covType)
{

    ProcessObservations();

    // Get the cross covariance between the evaluation points and the observations
    Eigen::MatrixXd crossCov;

    Eigen::MatrixXd outputMean(coDim, newLocs.cols());
    Eigen::MatrixXd outputCov;

    Eigen::Map<Eigen::VectorXd> postMean(outputMean.data(), coDim*newLocs.cols());

    if(observations.size()==0){
        outputMean = mean->Evaluate(newLocs);
        assert(outputMean.rows()==coDim);

    }else{
        // Get the cross covariance between the evaluation points and the observations
        crossCov = BuildCrossCov(newLocs);

        postMean = crossCov * sigmaTrainDiff;//covSolver.solve(trainDiff);
        outputMean += mean->Evaluate(newLocs);
    }

    // Only compute the prediction variances, not covariance
    if(covType == GaussianProcess::DiagonalCov){

        outputCov.resize(coDim, newLocs.cols());

        Eigen::MatrixXd priorCov(coDim, coDim);
        for(int i=0; i<newLocs.cols(); ++i){
            covKernel->FillCovariance(newLocs.col(i),priorCov);

            if(observations.size()>0){
                for(int d=0; d<coDim; ++d)
                    outputCov(d,i) = priorCov(d,d) - crossCov.row(i*coDim+d) * covSolver.solve(crossCov.row(i*coDim+d).transpose());
            }else{
                for(int d=0; d<coDim; ++d)
                    outputCov(d,i) = priorCov(d,d);
            }
        }

    // Predict marginal covariance at each point, but not the covariance between points
    }else if(covType==GaussianProcess::BlockCov){

        outputCov.resize(coDim, coDim*newLocs.cols());

        Eigen::MatrixXd priorCov(coDim, coDim);
        for(int i=0; i<newLocs.cols(); ++i){
            covKernel->FillCovariance(newLocs.col(i),priorCov);

            if(observations.size()>0){
                outputCov.block(0,coDim*i,coDim,coDim) = priorCov - crossCov.block(i*coDim,0,coDim,crossCov.cols()) * covSolver.solve(crossCov.block(i*coDim,0,coDim,crossCov.cols()).transpose());
            }else{
                outputCov.block(0,coDim*i,coDim,coDim) = priorCov;
            }
        }


    // Compute the full joint covariance of all predictions
    }else if(covType==GaussianProcess::FullCov){

        Eigen::MatrixXd priorCov(coDim*newLocs.cols(), coDim*newLocs.cols());
        covKernel->FillCovariance(newLocs,priorCov);

        // Solve for the posterior covariance
        if(observations.size()>0){
            outputCov = priorCov - crossCov * covSolver.solve(crossCov.transpose());
        }else{
            outputCov = priorCov;
        }
    }

    // Add a nugget to the covariance diagonal
    outputCov += nugget*Eigen::MatrixXd::Identity(outputCov.rows(), outputCov.cols());
    return std::make_pair(outputMean, outputCov);
};


Eigen::MatrixXd GaussianProcess::PredictMean(Eigen::MatrixXd const& newPts)
{
    if(observations.size()==0)
        return mean->Evaluate(newPts);

    ProcessObservations();

    // Construct the cross covariance
    Eigen::MatrixXd crossCov = BuildCrossCov(newPts);

    Eigen::MatrixXd meanMat(coDim, newPts.size());
    Eigen::Map<Eigen::VectorXd> meanVec(meanMat.data(), coDim*newPts.cols());

    meanVec = crossCov * sigmaTrainDiff;
    meanMat += mean->Evaluate(newPts);

    return meanMat;
}


Eigen::MatrixXd GaussianProcess::Sample(Eigen::MatrixXd const& newPts)
{
    ProcessObservations();

    Eigen::MatrixXd mean, covariance;
    std::tie(mean,covariance) = Predict(newPts, GaussianProcess::FullCov);

    Eigen::MatrixXd output(mean.rows(), mean.cols());
    Eigen::Map<Eigen::VectorXd> outVec(output.data(), output.rows()*output.cols());

    outVec = covariance.selfadjointView<Eigen::Lower>().llt().matrixL()*RandomGenerator::GetNormal(covariance.rows());

    output += mean;

    return output;
}

double GaussianProcess::LogLikelihood(Eigen::MatrixXd const& xs,
                                      Eigen::MatrixXd const& vals)
{

    ProcessObservations();

    Eigen::MatrixXd mean, covariance;
    std::tie(mean,covariance) = Predict(xs, GaussianProcess::FullCov);

    const Eigen::Map<const Eigen::VectorXd> valMap(vals.data(), vals.rows()*vals.cols());
    Eigen::Map<Eigen::VectorXd> muMap(mean.data(), mean.rows()*mean.cols());

    Eigen::VectorXd diff = valMap-muMap;
    auto solver = covariance.selfadjointView<Eigen::Lower>().llt();
    double logDet = 0;
    for(int i=0; i<solver.matrixL().rows(); ++i)
        logDet += 2.0*log(solver.matrixL()(i,i));

    double logDens = -0.5*diff.dot( solver.solve(diff) ) - 0.5*logDet - 0.5*valMap.size()*log(2.0*M_PI);
    return logDens;
}


double GaussianProcess::MarginalLogLikelihood()
{
    Eigen::VectorXd grad;
    return MarginalLogLikelihood(grad,false);
}

// Evaluates the log marginal likelihood needed when fitting hyperparameters
double  GaussianProcess::MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                               bool                        computeGrad)
{
  assert(false);
  //   ProcessObservations();
  //
  //   // Compute the log determinant of the covariance
  //   Eigen::MatrixXd L = covSolver.matrixL();
  //   L = (L*covSolver.vectorD().cwiseSqrt().asDiagonal()).eval();
  //   L = (covSolver.transpositionsP().transpose() * L).eval();
  //   double logDet = 0.0;
  //   for(int i=0; i<L.rows(); ++i)
	// logDet += 2.0*log(L(i,i));
  //
  //   // Make the mean prediction and compute the difference with observations
  //   double logLikely = -0.5 * trainDiff.transpose() * sigmaTrainDiff - 0.5*logDet - 0.5*observations.size()*log(2.0*pi);
  //
  //   if(computeGrad)
  //   {
	// // Compute the gradient of the log likelihood
	// const int numParams = covKernel->numParams;
	// grad.resize(numParams);
  //
	// for(int p=0; p<numParams; ++p)
	// {
  //           // Build the derivative matrix
  //           Eigen::MatrixXd derivMat(obsDim, obsDim);
  //           Eigen::MatrixXd tempDerivMat;
  //           int currRow=0;
  //           int currCol=0;
  //
  //           for(int j=0; j<observations.size(); ++j)
  //           {
  //               currRow = 0;
  //               for(int i=0; i<=j; ++i)
  //               {
  //                   tempDerivMat = covKernel->GetDerivative(observations.at(i)->loc, observations.at(j)->loc, p);
  //                   derivMat.block(currRow,currCol, observations.at(i)->H->rows(), observations.at(j)->H->rows()) = observations.at(i)->H->Apply( observations.at(j)->H->Apply(tempDerivMat).transpose() );
  //
  //                   currRow += observations.at(i)->H->rows();
  //               }
  //               currCol += observations.at(j)->H->rows();
  //           }
  //
	//     grad(p) = 0.5*(sigmaTrainDiff*sigmaTrainDiff.transpose()*derivMat - covSolver.solve(derivMat)).trace();
	// }
  //   }
  //
  //   return logLikely;
};
