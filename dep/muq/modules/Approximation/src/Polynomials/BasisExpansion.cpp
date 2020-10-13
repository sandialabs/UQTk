#include "MUQ/Approximation/Polynomials/BasisExpansion.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

using namespace muq::Approximation;
using namespace muq::Utilities;

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                               bool                                                    coeffInput) :
                                 BasisExpansion(basisCompsIn,
                                                MultiIndexFactory::CreateTotalOrder(basisCompsIn.size(),0))
{
}

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                               std::shared_ptr<muq::Utilities::MultiIndexSet>          multisIn,
                               bool                                                    coeffInput) :
                                 BasisExpansion(basisCompsIn,
                                                multisIn,
                                                Eigen::MatrixXd::Zero(1,multisIn->Size()))
{
};

BasisExpansion::BasisExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                               std::shared_ptr<muq::Utilities::MultiIndexSet>          multisIn,
                               Eigen::MatrixXd                                  const& coeffsIn,
                               bool                                                    coeffInput) :
                                 ModPiece(GetInputSizes(multisIn, coeffsIn, coeffInput), GetOutputSizes(multisIn, coeffsIn)),
                                 basisComps(basisCompsIn),
                                 multis(multisIn),
                                 coeffs(coeffsIn)
{
  assert(basisComps.size() == multis->GetMultiLength());
  assert(multis->Size() == coeffs.cols());
}

Eigen::VectorXi BasisExpansion::GetInputSizes(std::shared_ptr<muq::Utilities::MultiIndexSet> multisIn,
                                              Eigen::MatrixXd                         const& coeffsIn,
                                              bool                                           coeffInput)
{
  Eigen::VectorXi output;
  if(coeffInput){
    output.resize(2);
    output(0) = multisIn->GetMultiLength();
    output(1) = coeffsIn.cols()*coeffsIn.rows();
  }else{
    output.resize(1);
    output(0) = multisIn->GetMultiLength();
  }
  return output;
}

Eigen::VectorXi BasisExpansion::GetOutputSizes(std::shared_ptr<muq::Utilities::MultiIndexSet> multisIn,
                                               Eigen::MatrixXd                         const& coeffsIn)
{
  return coeffsIn.rows() * Eigen::VectorXi::Ones(1);
}


Eigen::VectorXd BasisExpansion::GetAllTerms(Eigen::VectorXd const& x) const{

  // Get the maximum orders
  Eigen::VectorXi maxOrders = multis->GetMaxOrders();

  // Evaluate each dimension up to the maximum order
  std::vector<std::vector<double>> uniEvals(basisComps.size());
  assert(uniEvals.size() == maxOrders.size());

  for(int i=0; i<uniEvals.size(); ++i){
    uniEvals.at(i).resize(maxOrders(i)+1);
    for(int j=0; j<=maxOrders(i); ++j){
      uniEvals.at(i).at(j) = basisComps.at(i)->BasisEvaluate(j, x(i));
    }
  }

  // Now that we have all the univariate terms evaluated, evaluate the expansion
  Eigen::VectorXd allTerms = Eigen::VectorXd::Ones(multis->Size());
  for(int i=0; i<multis->Size(); ++i){

    for(auto it = multis->at(i)->GetNzBegin(); it != multis->at(i)->GetNzEnd(); ++it)
      allTerms(i) *= uniEvals.at(it->first).at(it->second);
  }

  return allTerms;
}

Eigen::MatrixXd BasisExpansion::GetAllDerivs(Eigen::VectorXd const& x) const{

  // Get the maximum orders
  Eigen::VectorXi maxOrders = multis->GetMaxOrders();

  // Evaluate each dimension up to the maximum order
  std::vector<std::vector<double>> uniEvals(basisComps.size());
  std::vector<std::vector<double>> uniDerivs(basisComps.size());
  assert(uniEvals.size() == maxOrders.size());

  for(int i=0; i<uniEvals.size(); ++i){
    uniEvals.at(i).resize(maxOrders(i)+1);
    uniDerivs.at(i).resize(maxOrders(i)+1);

    for(int j=0; j<=maxOrders(i); ++j){
      uniEvals.at(i).at(j) = basisComps.at(i)->BasisEvaluate(j, x(i));
      uniDerivs.at(i).at(j) = basisComps.at(i)->DerivativeEvaluate(j, 1, x(i));
    }
  }

  // Now that we have all the univariate terms evaluated, evaluate the expansion
  Eigen::MatrixXd allDerivs = Eigen::MatrixXd::Ones(multis->Size(),x.size());
  for(int i=0; i<multis->Size(); ++i){

    // Loop over each dimension
    for(int j=0; j<x.size(); ++j){
      if(multis->at(i)->GetValue(j)==0){
        allDerivs(i,j) = 0;
      }else{
        for(auto it = multis->at(i)->GetNzBegin(); it != multis->at(i)->GetNzEnd(); ++it){

          if(it->first == j){
            allDerivs(i,j) *= uniDerivs.at(it->first).at(it->second);
          }else{
            allDerivs(i,j) *= uniEvals.at(it->first).at(it->second);
          }

        }
      }
    }
  }

  return allDerivs;
}

std::vector<Eigen::MatrixXd> BasisExpansion::GetHessians(Eigen::VectorXd const& x) const{

  // Get the maximum orders
  Eigen::VectorXi maxOrders = multis->GetMaxOrders();

  // Evaluate each dimension up to the maximum order
  std::vector<std::vector<double>> uniEvals(basisComps.size());
  std::vector<std::vector<double>> uniD1(basisComps.size()); // first derivatives
  std::vector<std::vector<double>> uniD2(basisComps.size()); // second derivatives

  assert(uniEvals.size() == maxOrders.size());

  for(int i=0; i<uniEvals.size(); ++i){
    uniEvals.at(i).resize(maxOrders(i)+1);
    uniD1.at(i).resize(maxOrders(i)+1);
    uniD2.at(i).resize(maxOrders(i)+1);

    for(int j=0; j<=maxOrders(i); ++j){
      uniEvals.at(i).at(j) = basisComps.at(i)->BasisEvaluate(j, x(i));
      uniD1.at(i).at(j) = basisComps.at(i)->DerivativeEvaluate(j, 1, x(i));
      uniD2.at(i).at(j) = basisComps.at(i)->DerivativeEvaluate(j, 2, x(i));
    }
  }

  std::vector<Eigen::MatrixXd> hessians(coeffs.rows(), Eigen::MatrixXd::Zero(x.size(),x.size()));

  // Loop over each term in the expansion
  for(int i=0; i<multis->Size(); ++i){

    // Loop over each dimension
    for(int j=0; j<x.size(); ++j){
      for(int k=0; k<=j; ++k){
        if((multis->at(i)->GetValue(j)!=0) && ((multis->at(i)->GetValue(k)!=0))){

          double tempVal = 1.0;

          for(auto it = multis->at(i)->GetNzBegin(); it != multis->at(i)->GetNzEnd(); ++it){

            if((j==k) && (it->first == j)){
              tempVal *= uniD2.at(it->first).at(it->second);
            }else if((it->first == j)||(it->first == k)){
              tempVal *= uniD1.at(it->first).at(it->second);
            }else{
              tempVal *= uniEvals.at(it->first).at(it->second);
            }
          }

          // Add the results into each of the hessians matrices
          for(int kk=0; kk<coeffs.rows(); ++kk)
            hessians.at(kk)(j,k) += coeffs(kk,i)*tempVal;
        }
      }
    }
  }

  // make sure all the hessians are symmetric
  for(int kk=0; kk<coeffs.rows(); ++kk)
    hessians.at(kk).triangularView<Eigen::Upper>() =  hessians.at(kk).triangularView<Eigen::Lower>().transpose();

  return hessians;

}

void BasisExpansion::ProcessCoeffs(Eigen::VectorXd const& newCoeffs)
{
  assert(newCoeffs.size() == coeffs.rows()*coeffs.cols());

  Eigen::Map<const Eigen::MatrixXd> coeffMap(newCoeffs.data(), coeffs.rows(), coeffs.cols());

  coeffs = coeffMap;
}


void BasisExpansion::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {

  Eigen::VectorXd const& x = inputs.at(0);
  if(inputs.size()>1)
    ProcessCoeffs(inputs.at(1));

  // Compute the output
  outputs.resize(1);
  outputs.at(0) = (coeffs*GetAllTerms(x)).eval();
}

void BasisExpansion::JacobianImpl(unsigned int const                           wrtOut,
                                  unsigned int const                           wrtIn,
                                  muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs)
{
  assert(wrtOut==0);

  Eigen::VectorXd const& x = inputs.at(0);
  if(inputs.size()>1)
    ProcessCoeffs(inputs.at(1));

  if(wrtIn==0){
    jacobian = Eigen::MatrixXd(coeffs*GetAllDerivs(x));
  }else if(wrtIn==1){
    jacobian = Eigen::MatrixXd(GetAllTerms(x).transpose().replicate(coeffs.rows(),1));
  }
}


Eigen::MatrixXd BasisExpansion::SecondDerivative(unsigned                       outputDim,
                                                 unsigned                       derivDim1,
                                                 unsigned                       derivDim2,
                                                 Eigen::VectorXd         const& x,
                                                 Eigen::MatrixXd         const& newCoeffs)
{
  SetCoeffs(newCoeffs);
  return SecondDerivative(outputDim, derivDim1, derivDim2, x);
}

Eigen::MatrixXd BasisExpansion::SecondDerivative(unsigned                       outputDim,
                                                 unsigned                       derivDim1,
                                                 unsigned                       derivDim2,
                                                 Eigen::VectorXd         const& x)
{
  if((derivDim1==0) && (derivDim2==1)){
    return GetAllDerivs(x).transpose();
  }else if((derivDim1==1) && (derivDim2==0)){
    return GetAllDerivs(x);
  }else if(derivDim1==0){
    return GetHessians(x).at(outputDim);
  }else{
    return Eigen::MatrixXd::Zero(coeffs.cols(), coeffs.cols());
  }
}


Eigen::MatrixXd BasisExpansion::GetCoeffs() const{
  return coeffs;
}

void BasisExpansion::SetCoeffs(Eigen::MatrixXd const& allCoeffs){
  assert(coeffs.rows()==allCoeffs.rows());
  assert(coeffs.cols()==allCoeffs.cols());
  coeffs = allCoeffs;
}

Eigen::MatrixXd BasisExpansion::BuildVandermonde(Eigen::MatrixXd const& evalPts) const
{
  Eigen::MatrixXd vand(evalPts.cols(), NumTerms());

  for(int i=0; i<evalPts.cols(); ++i)
    vand.row(i) = GetAllTerms(evalPts.col(i)).transpose();

  return vand;
}


Eigen::MatrixXd BasisExpansion::BuildDerivMatrix(Eigen::MatrixXd const& evalPts, int wrtDim) const
{
  Eigen::MatrixXd output(evalPts.cols(), NumTerms());

  for(int i=0; i<evalPts.cols(); ++i)
    output.row(i) = GetAllDerivs(evalPts.col(i)).col(wrtDim).transpose();

  return output;
}
