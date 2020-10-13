#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"


using namespace std;
using namespace muq::Utilities;

void muq::Utilities::MultiIndexFactory::RecursiveTotalOrderFill(unsigned int const                 maxOrder,
                                                                unsigned int const                 minOrder,
                                                                shared_ptr<MultiIndexSet>          output,
                                                                unsigned int const                 currDim,
                                                                Eigen::RowVectorXi                &base,
                                                                std::shared_ptr<MultiIndexLimiter> limiter)
{
  int currOrder = base.head(currDim+1).sum();
  const int length = base.size();

  if(currDim==length-1)
  {
    for(int i=max<int>(0,minOrder-currOrder); i<=maxOrder-currOrder; ++i)
    {
      base(length-1) = i;
      auto newTerm = make_shared<MultiIndex>(base);
      if(limiter->IsFeasible(newTerm))
        output->AddActive(newTerm);
    }
  }else{
    for(int i=0; i<=maxOrder-currOrder; ++i)
    {
      base.tail(length-currDim).setZero();
      base(currDim) = i;
      RecursiveTotalOrderFill(maxOrder,minOrder,output,currDim+1,base,limiter);
    }
  }
}

void muq::Utilities::MultiIndexFactory::RecursiveHyperbolicFill(const double                  maxNormPow,
                                                                shared_ptr<MultiIndexSet>     output,
                                                                unsigned int const            currDim,
                                                                Eigen::RowVectorXi           &base,
                                                                const double                  q,
                                                                shared_ptr<MultiIndexLimiter> limiter)
{
  double currNorm = 0;
  for(int i=0; i<currDim; ++i)
    currNorm += pow(static_cast<double>(base(i)),q);

  const int length = base.size();

  if(currDim==length-1)
  {
    double newNorm = currNorm;
    base(length-1) = 0;
    while(newNorm<maxNormPow)
    {
      auto newTerm = make_shared<MultiIndex>(base);
      if(limiter->IsFeasible(newTerm))
        output->AddActive(newTerm);
      base(length-1)++;
      newNorm = currNorm + pow(static_cast<double>(base(length-1)),q);
    }

  }else{
    double newNorm = currNorm;
    base.tail(length-currDim).setZero();
    while(newNorm<maxNormPow)
    {
      RecursiveHyperbolicFill(maxNormPow,output,currDim+1,base,q,limiter);
      base(currDim)++;
      newNorm = currNorm + pow(static_cast<double>(base(currDim)),q);
    }

  }
}


void muq::Utilities::MultiIndexFactory::RecursiveTensor(Eigen::RowVectorXi          const& orders,
                                                        shared_ptr<MultiIndexSet>          output,
                                                        unsigned int const                 currDim,
                                                        Eigen::RowVectorXi                &base,
                                                        std::shared_ptr<MultiIndexLimiter> limiter,
                                                        bool                               allInactive)
{
  const unsigned length = base.size();
  if(currDim==length-1)
  {
    int globalInd;
    shared_ptr<MultiIndex> newMulti;
    // add all the active indices
    for(int i=0; i<=orders(length-1); ++i)
    {
      base(length-1) = i;
      newMulti = make_shared<MultiIndex>(base);

      if(limiter->IsFeasible(newMulti)){
        if(allInactive){
          output->AddInactive(newMulti);
        }else{
          output->AddActive(newMulti);
        }
      }
    }

  }else{

    // add all the active nodes
    for(int i=0; i<=orders(currDim); ++i)
    {
      base.tail(length-currDim).setZero();
      base(currDim) = i;
      RecursiveTensor(orders,output,currDim+1,base,limiter, allInactive);
    }

    // add inactive neighboring nodes
    base.tail(length-currDim).setZero();
    base(currDim) = orders(currDim)+1;
    RecursiveTensor(orders,output,currDim+1,base,limiter, true);

  }
}


shared_ptr<MultiIndexSet> muq::Utilities::MultiIndexFactory::CreateTotalOrder(unsigned int const                 length,
                                                                              unsigned int const                 maxOrder,
                                                                              unsigned int const                 minOrder,
                                                                              std::shared_ptr<MultiIndexLimiter> limiter)
{
  assert(maxOrder>=minOrder);
  assert(minOrder>=0);
  assert(length>0);

  // create an empy multiindex set
  shared_ptr<MultiIndexSet> output = make_shared<MultiIndexSet>(length,limiter);

  // start with a vector of zeros
  Eigen::RowVectorXi base = Eigen::RowVectorXi::Zero(length);

  RecursiveTotalOrderFill(maxOrder,minOrder,output,0,base,limiter);

  return output;
}

std::vector<std::shared_ptr<MultiIndexSet>> muq::Utilities::MultiIndexFactory::CreateTriTotalOrder(unsigned int const                 length,
                                                                                                   unsigned int const                 maxOrder,
                                                                                                   unsigned int const                 minOrder,
                                                                                                   std::shared_ptr<MultiIndexLimiter> limiter)
{

  vector<shared_ptr<MultiIndexSet>> multis(length);
  for(int d=0; d<length; ++d){
    auto dimLimiter = make_shared<AndLimiter>(limiter,make_shared<DimensionLimiter>(0,d+1));
    multis.at(d) = MultiIndexFactory::CreateTotalOrder(length,maxOrder,minOrder,dimLimiter);
  }
  return multis;
}

shared_ptr<MultiIndexSet> muq::Utilities::MultiIndexFactory::CreateHyperbolic(unsigned int const                 length,
                                                                              unsigned int const                 maxOrder,
                                                                              const double                       q,
                                                                              std::shared_ptr<MultiIndexLimiter> limiter)
{
  assert(maxOrder>=0);
  assert(length>0);

  // create an empy multiindex set
  shared_ptr<MultiIndexSet> output = make_shared<MultiIndexSet>(length,limiter);

  // start with a vector of zeros
  Eigen::RowVectorXi base = Eigen::RowVectorXi::Zero(length);

  const double nugget = 1e-5;// needed for maximum power to be inclusive
  RecursiveHyperbolicFill(pow(static_cast<double>(maxOrder),q)+nugget,output,0,base,q,limiter);

  return output;
}

std::vector<std::shared_ptr<MultiIndexSet>> muq::Utilities::MultiIndexFactory::CreateTriHyperbolic(unsigned int const                 length,
                                                                                                   unsigned int const                 maxOrder,
                                                                                                   const double                       q,
                                                                                                   std::shared_ptr<MultiIndexLimiter> limiter)
{

  vector<shared_ptr<MultiIndexSet>> multis(length);
  for(int d=0; d<length; ++d){
    auto dimLimiter = make_shared<AndLimiter>(limiter,make_shared<DimensionLimiter>(0,d+1));
    multis.at(d) = MultiIndexFactory::CreateHyperbolic(length,maxOrder,q,dimLimiter);
  }
  return multis;
}

std::shared_ptr<MultiIndexSet> muq::Utilities::MultiIndexFactory::CreateFullTensor(unsigned int const length,
                                                                                   unsigned int const order,
                                                                                   std::shared_ptr<MultiIndexLimiter> limiter)
{
  return muq::Utilities::MultiIndexFactory::CreateFullTensor(order*Eigen::RowVectorXi::Ones(length), limiter);
}


std::shared_ptr<MultiIndexSet> muq::Utilities::MultiIndexFactory::CreateFullTensor(Eigen::RowVectorXi          const& orders,
                                                                                   std::shared_ptr<MultiIndexLimiter> limiter)
{
  assert(orders.minCoeff()>=0);
  assert(orders.size()>0);

  unsigned int length = orders.size();
  unsigned int numIndices = (orders.array()+2).prod(); // an overestimate of the total number of indices
  unsigned int numActiveIndices = (orders.array()+1).prod(); // the exact number of active indices

  // create an empy multiindex set
  shared_ptr<MultiIndexSet> output = make_shared<MultiIndexSet>(length,limiter);

  // start with a vector of zeros
  Eigen::RowVectorXi base = Eigen::RowVectorXi::Zero(length);

  RecursiveTensor(orders, output, 0, base, limiter, false);

  return output;
}


std::shared_ptr<MultiIndex> MultiIndexFactory::CreateSingleTerm(int totalDim, int nonzeroDim, int order)
{
  shared_ptr<MultiIndex> output = make_shared<MultiIndex>(totalDim);
  output->SetValue(nonzeroDim,order);
  return output;
}
