#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include <iostream>

using namespace muq::Approximation;
using namespace muq::Utilities;

FullTensorQuadrature::FullTensorQuadrature(unsigned int                       dim,
                                           std::shared_ptr<Quadrature> const& rule,
                                           unsigned int                       order) : FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>>(dim,rule),
                                                                                                            order*Eigen::RowVectorXi::Ones(dim)){};

FullTensorQuadrature::FullTensorQuadrature(unsigned int                       dim,
                                          std::shared_ptr<Quadrature> const& rule) : FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>>(dim,rule)){};


FullTensorQuadrature::FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>> const& rulesIn,
                                           Eigen::RowVectorXi                              orders) : Quadrature(rulesIn.size()),
                                                                                                       rules(rulesIn)
{
  for(int i=0; i<rules.size(); ++i)
    assert(rules.at(i)->Dim()==1);

  // If the quadrature orders are specified, call Compute
  if(orders.size()>0){
    assert(rulesIn.size()==orders.size());
    Compute(orders);
  }
}

void FullTensorQuadrature::Compute(unsigned int order)
{
  Compute(order*Eigen::RowVectorXi::Ones(dim));
}

void FullTensorQuadrature::Compute(Eigen::RowVectorXi const& orders) {

  assert(orders.size()==dim);

  Eigen::RowVectorXi tensorOrders(dim);

  // Compute the points and weights from the 1d rules
  std::vector<Eigen::MatrixXd> allPts(dim);
  std::vector<Eigen::VectorXd> allWts(dim);
  for(int i = 0; i<dim; ++i){
    rules.at(i)->Compute(orders(i));
    allPts.at(i) = rules.at(i)->Points();
    allWts.at(i) = rules.at(i)->Weights();

    tensorOrders(i) = allPts.at(i).cols()-1;
  }

  // First, compute all of the multiindices
  auto multis = MultiIndexFactory::CreateFullTensor(tensorOrders);

  // Now compute all the terms
  pts.resize(dim,multis->Size());
  wts = Eigen::VectorXd::Ones(multis->Size());

  for(int i=0; i<multis->Size(); ++i){
    auto& multi = multis->IndexToMulti(i);

    for(int d=0; d<dim; ++d){
      wts(i) *= allWts.at(d)(multi->GetValue(d));
      pts(d,i) = allPts.at(d)(0,multi->GetValue(d));
    }
  }

}

unsigned int FullTensorQuadrature::Exactness(unsigned int quadOrder) const {

  unsigned int maxOrder = 0;
  for(auto& rule : rules)
    maxOrder = std::max(maxOrder, rule->Exactness(quadOrder));
    
  return maxOrder;
}
