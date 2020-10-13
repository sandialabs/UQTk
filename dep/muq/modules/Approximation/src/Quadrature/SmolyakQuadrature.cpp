#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include "MUQ/Utilities/EigenUtilities.h"

#include <map>
#include <iostream>

using namespace muq::Approximation;
using namespace muq::Utilities;

SmolyakQuadrature::SmolyakQuadrature(unsigned int dim, std::shared_ptr<Quadrature> const& scalarRule) : SmolyakQuadrature(std::vector<std::shared_ptr<Quadrature>>(dim,scalarRule)){}


SmolyakQuadrature::SmolyakQuadrature(std::vector<std::shared_ptr<Quadrature>> const& scalarRulesIn) : Quadrature(scalarRulesIn.size()),
                                                                                               scalarRules(scalarRulesIn)
{}


void SmolyakQuadrature::Compute(unsigned int order)
{
  Compute(order*Eigen::RowVectorXi::Ones(dim));
}

void SmolyakQuadrature::Compute(Eigen::RowVectorXi const& orders)
{

  // Build the multiindex set.
  std::shared_ptr<MultiIndexSet> multis = BuildMultis(orders);

  Compute(multis);
}


std::shared_ptr<MultiIndexSet> SmolyakQuadrature::BuildMultis(Eigen::RowVectorXi const& orders) const
{
  // Use the minimum order in the vector to form a total order multiindex
  int minOrder = orders.minCoeff();
  assert(minOrder>=0);

  auto multis = MultiIndexFactory::CreateTotalOrder(dim,minOrder);

  // // Add other terms to get the right order
  // for(int i=0; i<dim; ++i){
  //   for(int p=minOrder+1; p<=orders(i); ++p)
  //   {
  //     auto newMulti = std::make_shared<MultiIndex>(dim);
  //     newMulti->SetValue(i,p);
  //     multis += newMulti;
  //   }
  // }

  return multis;
}


void SmolyakQuadrature::Compute(std::shared_ptr<MultiIndexSet> const& multis) {

  // Compute the weights caused by using a tensor product of quadrature rules directly
  Eigen::VectorXd smolyWts = ComputeWeights(multis);

  auto tensorQuad = std::make_shared<FullTensorQuadrature>(scalarRules);

  // A map holding the unique quadrature points (keys) and weights (values)
  std::map<Eigen::VectorXd, double, muq::Utilities::VectorLessThan<double>> quadParts;

  for(int i=0; i<multis->Size(); ++i)
  {

    if(std::abs(smolyWts(i))>5.0*std::numeric_limits<double>::epsilon()){
      Eigen::RowVectorXi multiVec = multis->IndexToMulti(i)->GetVector();

      tensorQuad->Compute(multiVec);

      auto& tensorPts = tensorQuad->Points();
      auto& tensorWts = tensorQuad->Weights();

      // Add all the tensor product points to the map of Smolyak points
      for(int ptInd = 0; ptInd<tensorPts.cols(); ++ptInd){

        auto iter = quadParts.find(tensorPts.col(ptInd));
        if(iter!=quadParts.end()){
          iter->second += smolyWts(i)*tensorWts(ptInd);
        }else{
          quadParts[tensorPts.col(ptInd)] = smolyWts(i)*tensorWts(ptInd);
        }
      }
    }
  }

  // Figure out how many nonzero weights we have
  unsigned int numNz = 0;
  double weightTol = 5.0*std::numeric_limits<double>::epsilon();
  for(auto& keyVal : quadParts)
    numNz += (std::abs(keyVal.second)>weightTol) ? 1.0 : 0.0;

  // Copy the unordered map, which has unique keys, into the Eigen matrix
  pts.resize(dim,numNz);
  wts.resize(numNz);
  unsigned int ind = 0;
  for(auto& part : quadParts){
    if(std::abs(part.second)>weightTol){
      pts.col(ind) = part.first;
      wts(ind) = part.second;
      ++ind;
    }
  }
}


void SmolyakQuadrature::UpdateWeights(unsigned int                          activeInd,
                                      std::shared_ptr<MultiIndexSet> const& multis,
                                      Eigen::VectorXd                     & multiWeights)
{
  unsigned int dim = multis->GetMultiLength();
  auto optMultis = MultiIndexFactory::CreateFullTensor(dim,1);

  auto& activeMulti = multis->IndexToMulti(activeInd);

  // This is the multiindex defining the Smolyak difference tensor product
  auto& k = multis->IndexToMulti(activeInd);

  for(unsigned int j=0; j<optMultis->Size(); ++j){

    auto newMulti = MultiIndex::Copy(activeMulti);
    double weightIncr = 1.0;

    // Loop over all the nonzero terms, which is where we have to decrement the order
    for(auto it = optMultis->IndexToMulti(j)->GetNzBegin(); it != optMultis->IndexToMulti(j)->GetNzEnd(); ++it) {
      if((it->second>0)&&(activeMulti->GetValue(it->first)>0)){
        weightIncr *= -1;
        newMulti->SetValue(it->first, activeMulti->GetValue(it->first)-1);
      }else if(k->GetValue(it->first)==0){
        weightIncr = 0.0;
      }
    }

    multiWeights(multis->MultiToIndex(newMulti)) += weightIncr;
  }
}

Eigen::VectorXd SmolyakQuadrature::ComputeWeights(std::shared_ptr<MultiIndexSet> const& multis)
{

  Eigen::VectorXd multiWeights = Eigen::VectorXd::Zero(multis->Size());

  for(unsigned int i = 0; i<multis->Size(); ++i) {
    UpdateWeights(i, multis, multiWeights);
  }

  return multiWeights;
}
