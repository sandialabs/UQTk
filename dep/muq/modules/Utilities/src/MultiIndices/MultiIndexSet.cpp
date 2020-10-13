#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include <algorithm>

using namespace muq::Utilities;
using namespace std;

std::shared_ptr<MultiIndexSet> muq::Utilities::operator+=( std::shared_ptr<MultiIndexSet> x,
                                                           std::shared_ptr<MultiIndexSet> y)
{
  (*x)+=(*y);
  return x;
}

std::shared_ptr<MultiIndexSet> muq::Utilities::operator+=( std::shared_ptr<MultiIndexSet> x,
                                                           std::shared_ptr<MultiIndex>    y)
{
  (*x)+=y;
  return x;
}

MultiIndexSet::MultiIndexSet(const unsigned dimIn,
                             shared_ptr<MultiIndexLimiter> limiterIn) : maxOrders(Eigen::VectorXi::Zero(dimIn)),
                                                                        dim(dimIn),
                                                                        limiter(limiterIn)
{
};

void MultiIndexSet::SetLimiter(std::shared_ptr<MultiIndexLimiter> const& limiterIn){

  // first, make sure no active terms in the set currently obey the new limiter.
  //  If a term is inactive, remove all edges tied to it
  for(int globalInd=0; globalInd<allMultis.size(); ++globalInd)
  {
    if(IsActive(globalInd)){
      if(!limiterIn->IsFeasible(allMultis.at(globalInd))){
        std::stringstream msg;
        msg << "Invalid limiter passed to MultiIndexSet::SetLimiter.  The active multi-index, ";
        msg << allMultis.at(globalInd)->GetVector() << ", is not valid with the new limiter.\n";
        throw std::invalid_argument(msg.str());
      }
    }else{

      if(!limiterIn->IsFeasible(allMultis.at(globalInd))){
        for(int inNode : inEdges[globalInd])
          outEdges[inNode].erase(globalInd);
        inEdges[globalInd].clear();
      }
    }
  }

  // copy the limiter
  limiter = limiterIn;
}

std::shared_ptr<MultiIndexSet> MultiIndexSet::CloneExisting(std::shared_ptr<MultiIndexSet> const& original)
{
  auto output = make_shared<MultiIndexSet>(original->dim, original->limiter);

  output->active2global  = original->active2global;
  output->outEdges       = original->outEdges;
  output->inEdges        = original->inEdges;
  output->maxOrders      = original->maxOrders;
  output->allMultis      = original->allMultis;

  return output;
}

// Eigen::MatrixXu MultiIndexSet::GetAllMultiIndices() const
// {
//
//   Eigen::MatrixXu output(active2local.size(), dim);
//   for(int i=0; i<active2local.size(); ++i){
//     auto multi = pool->at(local2global[active2local[i]]);
//     int multiDim = multi->GetDimension();
//     output.row(i).head(multiDim) = multi->GetMulti();
//   }
//
//   return output;
// }

int MultiIndexSet::MultiToIndex(std::shared_ptr<MultiIndex> const& input) const{

  auto localIter = multi2global.find(input);

  if(localIter!=multi2global.end()){
    return global2active[localIter->second];
  }else{
    return -1;
  }
}


int MultiIndexSet::AddMulti(std::shared_ptr<MultiIndex> const& newMulti)
{
  allMultis.push_back(newMulti);

  int globalInd = allMultis.size() - 1;
  multi2global[newMulti] = globalInd;

  global2active.push_back(-1);

  inEdges.push_back(std::set<int>());
  outEdges.push_back(std::set<int>());

  assert(allMultis.size() == global2active.size());

  AddForwardNeighbors(globalInd,false);
  AddBackwardNeighbors(globalInd, false);

  return globalInd;
}

///////////////////////////
//////// FINISHED UP HERE
///////////////////////////
int MultiIndexSet::AddActive(std::shared_ptr<MultiIndex> const& newNode)
{
  int globalInd = AddInactive(newNode);

  if(globalInd>=0){

    Activate(globalInd);
    return global2active[globalInd];

  }else{
    return -1;
  }
}



int MultiIndexSet::AddInactive(std::shared_ptr<MultiIndex> const& newNode)
{
  auto iter = multi2global.find(newNode);

  if(iter!=multi2global.end()){

    return iter->second;

  }else if(limiter->IsFeasible(newNode)){

    return AddMulti(newNode);

  }else{

    return -1;
  }
}

bool MultiIndexSet::IsActive(std::shared_ptr<MultiIndex> const& multiIndex) const
{
  auto iter = multi2global.find(multiIndex);

  if(iter!=multi2global.end()){
    return IsActive(iter->second);
  }else{
    return false;
  }
}

bool MultiIndexSet::IsActive(unsigned int globalIndex) const
{
  return global2active[globalIndex] >= 0;
}

bool MultiIndexSet::IsAdmissible(unsigned int globalIndex) const
{
  auto& multi = allMultis.at(globalIndex);

  if(!limiter->IsFeasible(multi))
    return false;

  if(IsActive(globalIndex))
    return true;

  // count the number of input edges that are coming from active indices
  int numAdmiss = 0;
  for(int inNode : inEdges.at(globalIndex)){
    if(IsActive(inNode))
      numAdmiss++;
  }

  if(numAdmiss==multi->nzInds.size()){
    return true;
  }else{
    return false;
  }
}

bool MultiIndexSet::IsAdmissible(std::shared_ptr<MultiIndex> const& multiIndex) const
{
  auto iter = multi2global.find(multiIndex);
  if(iter==multi2global.end()){
    return false;
  }else{
    return IsAdmissible(iter->second);
  }
}


bool MultiIndexSet::IsExpandable(unsigned int activeIndex) const
{
  // an index is expandable when at least one forward neighbor is admissible but not active (i.e. outedge)

  // loop through the outgoing edges for this node
  for(int nextInd : outEdges[active2global.at(activeIndex)]){
    if(!IsActive(nextInd)&&IsAdmissible(nextInd))
      return true;
  }
  return false;
}

void MultiIndexSet::Activate(int globalIndex)
{

  // the index is already in the global set, if the value is non-negative, it is also active and we don't need to do anything
  if(global2active.at(globalIndex)<0)
  {
    auto& newNode = allMultis.at(globalIndex);

    // now add the index to the active set
    active2global.push_back(globalIndex);

    int newActiveInd = active2global.size()-1;

    global2active.at(globalIndex) = newActiveInd;

    // update the maximum order
    for(auto pair : newNode->nzInds)
      maxOrders(pair.first) = std::max<unsigned>(maxOrders(pair.first),pair.second);

    AddForwardNeighbors(globalIndex,true);
    AddBackwardNeighbors(globalIndex,true);
  }
}

void MultiIndexSet::Activate(std::shared_ptr<MultiIndex> const& multiIndex)
{
  auto iter = multi2global.find(multiIndex);

  assert(iter!=multi2global.end());
  assert(IsAdmissible(iter->second));

  Activate(iter->second);
}

void MultiIndexSet::AddForwardNeighbors(unsigned int globalIndex, bool addInactive)
{

  Eigen::RowVectorXi base = allMultis.at(globalIndex)->GetVector();

  shared_ptr<MultiIndex> newNode;
  for(unsigned int i=0; i<base.size(); ++i)
  {
    base(i)++;

    newNode = make_shared<MultiIndex>(base);

    if(limiter->IsFeasible(newNode)){

      auto iter = multi2global.find(newNode);

      if(iter!=multi2global.end()){
        inEdges.at(iter->second).insert(globalIndex);
        outEdges.at(globalIndex).insert(iter->second);
      }else if(addInactive){
        AddInactive(newNode);
      }
    }
    base(i)--;
  }
}

std::vector<shared_ptr<MultiIndex>>  MultiIndexSet::GetAdmissibleForwardNeighbors(unsigned int activeIndex)
{
  unsigned int globalInd = active2global.at(activeIndex);

  vector<shared_ptr<MultiIndex>> output;
  for( auto neighbor : outEdges[globalInd])
  {
    if(IsAdmissible(neighbor))
      output.push_back(allMultis.at(neighbor));
  }

  return output;
}

std::vector<unsigned int> MultiIndexSet::GetFrontier() const {

  std::vector<unsigned int> frontierInds;

  for(unsigned int activeInd = 0; activeInd<active2global.size(); ++activeInd) {
    if(IsExpandable(activeInd))
      frontierInds.push_back(activeInd);
  }

  return frontierInds;
}

std::vector<unsigned int> MultiIndexSet::GetStrictFrontier() const
{
  std::vector<unsigned int> frontierInds;

  for(unsigned int activeInd = 0; activeInd<active2global.size(); ++activeInd) {
    // loop over all the forward neighbors
    unsigned int globalInd = active2global.at(activeInd);

    bool isStrict = true;
    for( auto neighbor : outEdges[globalInd])
      isStrict = isStrict && (!IsActive(neighbor));

    if(isStrict)
      frontierInds.push_back(activeInd);
  }

  return frontierInds;
}

std::vector<unsigned int> MultiIndexSet::GetBackwardNeighbors(unsigned int activeIndex) const
{
  unsigned int globalInd = active2global.at(activeIndex);

  std::vector<unsigned int> output;
  for(auto neighbor : inEdges[globalInd])
    output.push_back(global2active.at(neighbor));

  return output;
}

std::vector<unsigned int> MultiIndexSet::GetBackwardNeighbors(std::shared_ptr<MultiIndex> const& multiIndex) const
{
  auto iter = multi2global.find(multiIndex);

  assert(iter!=multi2global.end());

  unsigned int globalInd = iter->second;
  std::vector<unsigned int> output;
  for(auto neighbor : inEdges[globalInd])
    output.push_back(global2active.at(neighbor));

  return output;
}


unsigned int MultiIndexSet::NumActiveForward(unsigned int activeInd) const
{
  unsigned int globalInd = active2global.at(activeInd);

  unsigned int numActive = 0;
  for( auto neighbor : outEdges[globalInd])
  {
    if(IsActive(neighbor))
      numActive++;
  }
  return numActive;
}

unsigned int MultiIndexSet::NumForward(unsigned int activeInd) const
{
  unsigned int globalInd = active2global.at(activeInd);
  return outEdges[globalInd].size();
}

void MultiIndexSet::AddBackwardNeighbors(unsigned int globalIndex, bool addInactive)
{
  Eigen::RowVectorXi base = allMultis.at(globalIndex)->GetVector();

  shared_ptr<MultiIndex> newNode;
  for(unsigned int i=0; i<base.size(); ++i)
  {
    if(base(i)>0){

      base(i)--;

      newNode = make_shared<MultiIndex>(base);

      if(limiter->IsFeasible(newNode)){
        auto iter = multi2global.find(newNode);

        if(iter!=multi2global.end()){
          outEdges.at(iter->second).insert(globalIndex);
          inEdges.at(globalIndex).insert(iter->second);
        }else if(addInactive){
          AddInactive(newNode);
        }
      }
      base(i)++;
    }
  }
}

std::vector<unsigned> MultiIndexSet::Expand(unsigned int activeIndex)
{
  if(activeIndex >= active2global.size()){
    std::stringstream msg;
    msg << "Invalid index passed to MultiIndexSet::Expand.  A value of " << activeIndex << " was passed to the function, but only " << active2global.size() << " active components exist in the set.\n";
    throw std::out_of_range(msg.str());
  }

  vector<unsigned> newIndices;
  unsigned globalIndex = active2global.at(activeIndex);

  // loop through the forward neighbors of this index
  std::set<int> tempSet = outEdges.at(globalIndex);
  for(int neighbor : tempSet)
  {
    if(IsAdmissible(neighbor)&&(!IsActive(neighbor))){
      Activate(neighbor);
      newIndices.push_back(global2active.at(neighbor));
    }
  }

  // return the vector of newly activated indices
  return newIndices;
}

std::vector<unsigned> MultiIndexSet::ForciblyExpand(unsigned int const activeIndex)
{
  assert(activeIndex<active2global.size());

  vector<unsigned> newIndices;
  unsigned globalIndex = active2global.at(activeIndex);

  // loop through the forward neighbors of this index
  std::set<int>& tempSet = outEdges.at(globalIndex);
  for(int neighbor : tempSet)
    ForciblyActivate(neighbor,newIndices);

  // return the vector of newly activated indices
  return newIndices;

}

void MultiIndexSet::ForciblyActivate(int globalIndex, std::vector<unsigned> &newIndices){


  if(!IsActive(globalIndex)){

    // make the node active and add inactive neighbors if necessary, this also updates the edges and enables the loop below
    Activate(globalIndex);
    newIndices.push_back(global2active.at(globalIndex));

    // now, fill in all of the previous neighbors
    std::set<int>& tempSet = inEdges.at(globalIndex);
    for(int ind : tempSet)
      ForciblyActivate(ind,newIndices);

  }
}

std::vector<unsigned> MultiIndexSet::ForciblyActivate(std::shared_ptr<MultiIndex> const& multiIndex){

  assert(limiter->IsFeasible(multiIndex));

  auto iter = multi2global.find(multiIndex);
  vector<unsigned int> newIndices;

  // if we found the multiindex and it is active, there is nothing to do
  if(iter!=multi2global.end()){
    ForciblyActivate(iter->second,newIndices);
  }else{
    // Add the new index as an active node
    int newGlobalInd = AddInactive(multiIndex);
    ForciblyActivate(newGlobalInd,newIndices);
  }

  return newIndices;
}

MultiIndexSet& MultiIndexSet::operator+=(const MultiIndexSet& rhs)
{
  Union(rhs);
  return *this;
}

int MultiIndexSet::Union(const MultiIndexSet& rhs)
{
  int oldTerms = Size();

  for(int i = 0; i < rhs.allMultis.size(); ++i) {

    auto newMulti = rhs.allMultis.at(i);
    if(limiter->IsFeasible(newMulti)){
      if(rhs.global2active[i]<0){
        AddInactive(newMulti);
      }else{
        AddActive(newMulti);
      }
    }
  }

  return Size() - oldTerms;
}

MultiIndexSet& MultiIndexSet::operator+=(std::shared_ptr<MultiIndex> const& rhs)
{
  AddActive(rhs);
  return *this;
}
