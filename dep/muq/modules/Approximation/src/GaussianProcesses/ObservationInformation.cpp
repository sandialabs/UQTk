#include "MUQ/Approximation/GaussianProcesses/ObservationInformation.h"


using namespace muq::Utilities;
using namespace muq::Approximation;


void ObservationInformation::FillSelfCov(std::shared_ptr<KernelBase> kernel,
                                         Eigen::Ref<Eigen::MatrixXd> covBlock)
{
  covBlock = H->Apply( H->Apply( BuildBaseCovariance(kernel) ).transpose() );
}

void ObservationInformation::FillCrossCov(Eigen::Ref<const Eigen::VectorXd> const& otherLoc,
                                          std::shared_ptr<KernelBase>              kernel,
                                          Eigen::Ref<Eigen::MatrixXd>              covBlock)
{
  covBlock = H->Apply( BuildBaseCovariance(otherLoc, kernel) );
}

void ObservationInformation::FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                          std::shared_ptr<KernelBase>              kernel,
                                          Eigen::Ref<Eigen::MatrixXd>             covBlock)
{
    covBlock = otherObs->H->Apply( H->Apply( BuildBaseCovariance(otherObs, kernel) ).transpose() ).transpose();
};

Eigen::MatrixXd ObservationInformation::BuildBaseCovariance(Eigen::Ref<const Eigen::VectorXd> const& otherLoc,
                                                            std::shared_ptr<KernelBase>              kernel)
{
  return kernel->BuildCovariance(loc, otherLoc);
}

Eigen::MatrixXd ObservationInformation::BuildBaseCovariance(std::shared_ptr<KernelBase> kernel)
{
  return kernel->BuildCovariance(loc,loc);
}

Eigen::MatrixXd ObservationInformation::BuildBaseCovariance(std::shared_ptr<ObservationInformation> otherObs,
                                                            std::shared_ptr<KernelBase>              kernel)
{
  // Check to see if the other observation is a derivative observation
  auto derivObs = std::dynamic_pointer_cast<DerivativeObservation>(otherObs);
  if(derivObs){
    return derivObs->BuildBaseCovariance(shared_from_this(), kernel).transpose();
  }else{
    return kernel->BuildCovariance(loc, otherObs->loc);
  }
}



Eigen::MatrixXd DerivativeObservation::BuildBaseCovariance(Eigen::Ref<const Eigen::VectorXd> const& otherLoc,
                                                           std::shared_ptr<KernelBase>              kernel)
{
  Eigen::MatrixXd output(derivCoords.size() * kernel->coDim, kernel->coDim);
  for(int i=0; i<derivCoords.size(); ++i)
    output.block(i*derivCoords.size(),0,kernel->coDim,kernel->coDim) = kernel->GetPosDerivative(otherLoc, loc, derivCoords.at(i));

  return output;
}

Eigen::MatrixXd DerivativeObservation::BuildBaseCovariance(std::shared_ptr<KernelBase> kernel)
{
  Eigen::MatrixXd output(derivCoords.size() * kernel->coDim, derivCoords.size() * kernel->coDim);
  for(int i=0; i<derivCoords.size(); ++i){
    for(int j=0; j<derivCoords.size(); ++j){

      std::vector<int> allWrts = derivCoords.at(i);
      allWrts.insert(allWrts.end(), derivCoords.at(j).begin(), derivCoords.at(j).end());
      output.block(i*derivCoords.size(),j*derivCoords.size(),kernel->coDim,kernel->coDim) = kernel->GetPosDerivative(loc, loc, allWrts);
    }
  }

  return output;

}

Eigen::MatrixXd DerivativeObservation::BuildBaseCovariance(std::shared_ptr<ObservationInformation> otherObs,
                                                           std::shared_ptr<KernelBase>              kernel)
{

  // Check to see if the other observation is a derivative observation
  auto derivObs = std::dynamic_pointer_cast<DerivativeObservation>(otherObs);
  if(derivObs){

    Eigen::MatrixXd output(derivCoords.size() * kernel->coDim, derivObs->derivCoords.size() * kernel->coDim);

    for(int i=0; i<derivCoords.size(); ++i){
      for(int j=0; j<derivObs->derivCoords.size(); ++j){

        std::vector<int> allWrts = derivCoords.at(i);
        allWrts.insert(allWrts.end(), derivObs->derivCoords.at(j).begin(), derivObs->derivCoords.at(j).end());
        output.block(i*derivCoords.size(),j*derivObs->derivCoords.size(),kernel->coDim,kernel->coDim) = kernel->GetPosDerivative(loc, otherObs->loc, allWrts);
      }
    }

    return output;

  }else{
    return BuildBaseCovariance(otherObs->loc, kernel);
  }
}
