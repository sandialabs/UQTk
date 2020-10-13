#include "MUQ/Approximation/PolynomialChaos/SmolyakEstimator.h"

#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"
#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"

#include <chrono>

using namespace muq::Modeling;
using namespace muq::Utilities;
using namespace muq::Approximation;

template<typename EstimateType>
SmolyakEstimator<EstimateType>::SmolyakEstimator(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn)
  : model(modelIn),
    termMultis(std::make_shared<MultiIndexSet>(modelIn->inputSizes(0))),
    pointCache(modelIn->inputSizes(0))
{
  // make sure the model has a single input and output
  assert(model->inputSizes.size()==1);
  assert(model->outputSizes.size()==1);
}

template<typename EstimateType>
EstimateType SmolyakEstimator<EstimateType>::ComputeWeightedSum(Eigen::VectorXd const& weights) const
{
  assert(weights.size()<=terms.size());

  unsigned int firstNzInd = 0;
  for(unsigned int i=0; i<weights.size(); ++i){
    if(std::abs(weights(i))>nzTol) {
      firstNzInd = i;
      break;
    }
  }

  assert(std::abs(weights(firstNzInd))>nzTol);

  // compute the number of nonzero terms
  EstimateType res = AddEstimates(0.0, terms.at(firstNzInd).val, weights(firstNzInd), terms.at(firstNzInd).val);
  for(unsigned int i=firstNzInd+1; i<weights.size(); ++i){
    if(std::abs(weights(i))>nzTol)
      res = AddEstimates(1.0, res, weights(i), terms.at(i).val);
  }

  return res;
}


template<typename EstimateType>
EstimateType SmolyakEstimator<EstimateType>::ComputeWeightedSum() const
{
  unsigned int firstNzInd = 0;
  const double weightTol = 10.0*std::numeric_limits<double>::epsilon();
  for(unsigned int i=0; i<terms.size(); ++i){
    if(std::abs(terms.at(i).weight)>weightTol) {
      firstNzInd = i;
      break;
    }
  }

  // compute the number of nonzero terms
  EstimateType res = AddEstimates(0.0, terms.at(firstNzInd).val, terms.at(firstNzInd).weight, terms.at(firstNzInd).val);
  for(unsigned int i=firstNzInd+1; i<terms.size(); ++i){
    if(std::abs(terms.at(i).weight)>weightTol){
      res = AddEstimates(1.0, res, terms.at(i).weight, terms.at(i).val);
    }
  }

  return res;
}

template<typename EstimateType>
EstimateType SmolyakEstimator<EstimateType>::Adapt(boost::property_tree::ptree options)
{
  auto t_start = std::chrono::high_resolution_clock::now();
  double runtime = 0.0; // <- Time since the beginning of adapting (in seconds)

  UpdateErrors();

  timeTol = options.get("MaximumAdaptTime", std::numeric_limits<double>::infinity());
  errorTol = options.get("ErrorTol", 0.0);
  maxNumEvals = options.get("MaximumEvals",std::numeric_limits<unsigned int>::max());

  // If none of the stopping criteria have been set, throw an error
  if((timeTol==std::numeric_limits<double>::infinity())&&(errorTol==0)&&(maxNumEvals==std::numeric_limits<unsigned int>::max())){
    std::stringstream msg;
    msg << "Error in SmolyakEstimator::Adapt.  No stopping criteria were set.  ";
    msg << "The options property_tree must specify at least one of the following ";
    msg << "variables, \"MaximumAdaptTime\", \"ErrorTol\", or \"MaximumEvals\".";
    throw std::runtime_error(msg.str());
  }

  // Adapt until we've reached a stopping criteria
  while((globalError>errorTol)&&(runtime<timeTol)&&(numEvals<maxNumEvals)) {
    Refine();
    UpdateErrors();

    auto t_curr =std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration<double, std::milli>(t_curr-t_start).count()/1e3;

    // keep track of the convergence diagnostics
    errorHistory.push_back(globalError);
    timeHistory.push_back(runtime+timeHistory.at(0)); // Adding the zero is necessary because the first time corresponds to terms added before adaptation
  }


  return ComputeWeightedSum();
}

template<typename EstimateType>
bool SmolyakEstimator<EstimateType>::Refine()
{
  // Of all the terms on the frontier, find the one with the largest error indicator and expand it
  int expandInd = -1;
  double maxError = 0.0;
  for(unsigned int termInd : termMultis->GetFrontier()){

    // Make sure all the backward neighbors are old
    if(!terms.at(termInd).isOld){
      if(terms.at(termInd).localError > maxError){
        expandInd = termInd;
        maxError = terms.at(termInd).localError;
      }
    }
  }

  // If no admissible terms were found, just return
  if(expandInd==-1)
    return false;

  terms.at(expandInd).isOld = true;

  // Add the new terms to the set
  std::vector<std::shared_ptr<MultiIndex>> neighs = termMultis->GetAdmissibleForwardNeighbors(expandInd);

  std::vector<std::shared_ptr<MultiIndex>> termsToAdd;
  for(auto& neigh : neighs){

    // Make sure the forward neighbor is not already in the set and has old backward neighbors
    if(!termMultis->IsActive(neigh)){

      bool shouldAdd = true;

      // Look through all of the backward neighbors. They should all be "old" to add this term.
      std::vector<unsigned int> backNeighInds = termMultis->GetBackwardNeighbors(neigh);

      for(auto backInd : backNeighInds)
        shouldAdd = shouldAdd && terms.at(backInd).isOld;

      if(shouldAdd)
        termsToAdd.push_back(neigh);
    }


  }

  AddTerms(termsToAdd);

  return true;
}

template<typename EstimateType>
void SmolyakEstimator<EstimateType>::AddTerms(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet)
{
  std::vector<std::shared_ptr<MultiIndex>> multiVec(fixedSet->Size());
  for(unsigned int i=0; i<fixedSet->Size(); ++i)
    multiVec.at(i) = fixedSet->IndexToMulti(i);

  AddTerms(multiVec);
}

template<typename EstimateType>
std::vector<std::vector<Eigen::VectorXd>> SmolyakEstimator<EstimateType>::PointHistory() const
{
  std::vector<std::vector<Eigen::VectorXd>> output(pointHistory.size());
  for(int adaptInd=0; adaptInd<pointHistory.size(); ++adaptInd)
  {
    output.at(adaptInd).reserve(pointHistory.at(adaptInd).size());
    for(auto& ptInd : pointHistory.at(adaptInd))
      output.at(adaptInd).push_back( GetFromCache(ptInd) );
  }
  return output;
}


template<typename EstimateType>
void SmolyakEstimator<EstimateType>::AddTerms(std::vector<std::shared_ptr<MultiIndex>> const& fixedSet)
{
  /* For each term in the fixed set we want to construct a tensor product approximation.
     The i^th tensor product approximation will require model evaluations at N_i
     points.   The allNewPts vector will store all of the points that need to be evaluated
     and stored in our cache.  The evalInds vector will keep track of which points
     were needed for each tensor product approximation by storing the index of
     the cached points.
  */


  termHistory.push_back(fixedSet);

  for(unsigned int i=0; i<fixedSet.size(); ++i) {
    unsigned int termInd = terms.size();
    terms.push_back(SmolyTerm());
    termMultis->AddActive(fixedSet.at(i));

    assert(terms.size()==termMultis->Size());


    std::vector<Eigen::VectorXd> pts = OneTermPoints(fixedSet.at(i));

    // Figure out if we've already evaluated a point, or if we need to
    terms.at(termInd).evalInds.resize(pts.size());
    for(unsigned int k=0; k<pts.size(); ++k){

      // Check if the point is already in the cache
      int cacheId = InCache(pts.at(k));
      if(cacheId<0){
        int cacheInd = AddToCache(pts.at(k));
        terms.at(termInd).evalInds.at(k) = cacheInd;
        evalCache.push_back( Eigen::VectorXd() );

      }else{
        terms.at(termInd).evalInds.at(k) = cacheId;
      }
    }
  }

  // Update the weights
  Eigen::VectorXd newWeights = SmolyakQuadrature::ComputeWeights(termMultis);
  for(unsigned int j=0; j<newWeights.size(); ++j)
    terms.at(j).weight = newWeights(j);

  // Compute differential weights for the frontier terms
  for(unsigned int termInd : termMultis->GetFrontier()){
    terms.at(termInd).diffWeights = Eigen::VectorXd::Zero(terms.size());
    SmolyakQuadrature::UpdateWeights(termInd, termMultis, terms.at(termInd).diffWeights);
  }

  ////////////////////////////////////////////////////
  // Figure out which terms we need to actually compute.

  // All terms with nonzero Smolyak weights
  for(unsigned int termInd=0; termInd<terms.size(); ++termInd){
    if(std::abs(terms.at(termInd).weight) > nzTol)
      terms.at(termInd).isNeeded = true;
  }

  // All terms needed to compute error estimates on the "frontier" terms
  for(unsigned int termInd : termMultis->GetFrontier()){
    for(unsigned int wInd=0; wInd<terms.at(termInd).diffWeights.size(); ++wInd){
      if(std::abs(terms.at(termInd).diffWeights(wInd))>nzTol)
        terms.at(wInd).isNeeded = true;
    }
  }

  ////////////////////////////////////////////////////////
  // Figure out what points we need to evaluate to build the terms we need
  std::set<unsigned int> ptsToEval;
  for(unsigned int termInd=0; termInd<terms.size(); ++termInd) {

    // If the smolyak weight is zero or the term has already been computed, don't bother computing any needed points
    if((terms.at(termInd).isNeeded) && (!terms.at(termInd).isComputed)) {
      for(unsigned int i=0; i<terms.at(termInd).evalInds.size(); ++i) {
        unsigned int ptInd = terms.at(termInd).evalInds.at(i);
        // If the size of the output is zero, we haven't evaluated the model at this point yet
        if(evalCache.at(ptInd).size()==0)
          ptsToEval.insert(ptInd);
      }
    }
  }

  ////////////////////////////////////////////////////////
  // Evaluate all the points we need to -- could easily be parallelized
  EvaluatePoints(ptsToEval);

  pointHistory.push_back(ptsToEval);
  evalHistory.push_back(numEvals);

  ////////////////////////////////////////////////////////
  // Now, compute any new tensor product estimates that are necessary
  for(unsigned int i=0; i<termMultis->Size(); ++i) {

    // If we haven't built this estimate before, do it now
    if((terms.at(i).isNeeded) && (!terms.at(i).isComputed)) {

      // Copy references of the model output to a vector
      std::vector<std::reference_wrapper<const Eigen::VectorXd>> evals;
      for(unsigned int ptInd=0; ptInd<terms.at(i).evalInds.size(); ++ptInd){
        evals.push_back( evalCache.at(terms.at(i).evalInds.at(ptInd) ) );
      }
      terms.at(i).val = ComputeOneTerm(termMultis->IndexToMulti(i), evals);
      terms.at(i).isComputed = true;
    }

  }

}

template<typename EstimateType>
void SmolyakEstimator<EstimateType>::EvaluatePoints(std::set<unsigned int> const& ptsToEval)
{
  for(auto& ptInd : ptsToEval){
    evalCache.at(ptInd) = model->Evaluate( GetFromCache(ptInd) ).at(0);
    numEvals++;
  }
}

template<typename EstimateType>
void SmolyakEstimator<EstimateType>::UpdateErrors()
{
  globalError = 0.0;
  // Consider terms that have non-active forward neighbors, which means it can be expandable
  std::vector<unsigned int> frontierTerms = termMultis->GetFrontier();

  // For each frontier term, compute a local error estimate
  for(unsigned int termInd : frontierTerms){
    if(!terms.at(termInd).isOld){
      terms.at(termInd).localError = ComputeMagnitude( ComputeWeightedSum( terms.at(termInd).diffWeights ));
      globalError += terms.at(termInd).localError;
    }
  }
}

template<typename EstimateType>
EstimateType SmolyakEstimator<EstimateType>::Compute(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet,
                                                     boost::property_tree::ptree                           options)
{
  auto t_start = std::chrono::high_resolution_clock::now();
  double runtime = 0.0; // <- Time since the beginning of adapting (in seconds)

  Reset();

  AddTerms(fixedSet);


  // Make all terms not on the strict frontier old
  for(unsigned i=0;i<terms.size();++i)
    terms.at(i).isOld = true;
  for(unsigned int termInd : termMultis->GetStrictFrontier())
    terms.at(termInd).isOld = false;

  UpdateErrors();

  // Keep track of how long that took us
  auto t_curr =std::chrono::high_resolution_clock::now();
  runtime = std::chrono::duration<double, std::milli>(t_curr-t_start).count()/1e3;
  timeHistory.push_back(runtime);
  errorHistory.push_back(globalError);

  bool shouldAdapt = options.get("ShouldAdapt",false);
  if(shouldAdapt)
    return Adapt(options);

  // We've done all the work, just return a weighted sum of the tensor product approximations
  return ComputeWeightedSum();
}

template<typename EstimateType>
void SmolyakEstimator<EstimateType>::Reset()
{
  termMultis = std::make_shared<MultiIndexSet>(model->inputSizes(0));
  terms.clear();
  globalError = std::numeric_limits<double>::infinity();

  errorHistory.clear();
  evalHistory.clear();
  timeHistory.clear();
  termHistory.clear();
  pointHistory.clear();
  numEvals = 0;
}


// template<typename EstimateType>
// void SmolyakEstimator<EstimateType>::UpdateSmolyCoeffs(unsigned int const index)
// {
//   // Get backward neighbors
//   std::vector<unsigned int> neighInds = termMultis->GetBackwardNeighbors(index);
//
//   Eigen::RowVectorXi baseMulti = termMultis->IndexToMulti(index)->GetVector();
//   smolyWeights.at(index) += 1;
//
//   for(unsigned int neighInd : neighInds) {
//     int parity = (baseMulti - termMultis->IndexToMulti(neighInd)->GetVector()).sum();
//     smolyWeights.at(neighInd) += ( (parity % 2 == 0) ? 1 : -1);
//   }
// }

template<typename EstimateType>
int SmolyakEstimator<EstimateType>::AddToCache(Eigen::VectorXd const& newPt)
{
  int cacheId = InCache(newPt);
  if(cacheId<0){
    pointCache.add(newPt);
    return CacheSize()-1;
  }else{
    return cacheId;
  }
}

template<typename EstimateType>
int SmolyakEstimator<EstimateType>::SmolyakEstimator::InCache(Eigen::VectorXd const& input) const
{
  if(CacheSize()>0){
    std::vector<size_t> indices;
    std::vector<double> squaredDists;
    std::tie(indices, squaredDists) = pointCache.query(input, 1);

    if(squaredDists.at(0)<cacheTol){
      return indices.at(0);
    }else{
      return -1;
    }

  }else{
    return -1;
  }
}

namespace muq{
namespace Approximation{
  template class SmolyakEstimator<Eigen::VectorXd>;
  template class SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>;
}
}
