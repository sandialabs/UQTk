#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Modeling;

FlannCache::FlannCache(std::shared_ptr<ModPiece> function) : ModPiece(function->inputSizes, function->outputSizes), // can only have one input and output
							     function(function),
							     kdTree(std::make_shared<DynamicKDTreeAdaptor<>>(function->inputSizes(0))) {

  // the target function can only have one input/output
  assert(function->numInputs==1);
  assert(function->numOutputs==1);
	centroid = Eigen::VectorXd::Zero(inputSizes(0));
}

FlannCache::~FlannCache() {}

void FlannCache::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
    int cacheId = InCache(inputs.at(0));
    outputs.resize(1);
    if(cacheId < 0){
      Add(inputs.at(0));
      outputs.at(0) = outputCache.at(outputCache.size()-1);
    }else{
      outputs.at(0) = outputCache.at(cacheId);
    }
}

int FlannCache::InCache(Eigen::VectorXd const& input) const {
  if( Size()>0 ) { // if there are points in the cache
    std::vector<size_t> indices;
    std::vector<double> squaredDists;
    std::tie(indices, squaredDists) = kdTree->query(input, 1);

    if(squaredDists.at(0)<std::numeric_limits<double>::epsilon()){
      return indices.at(0);
    }
  }

  // the cache is either empty or none of the points in a small radius are exactly the point we care about
  return -1;
}

Eigen::VectorXd FlannCache::Add(Eigen::VectorXd const& newPt) {
  // evaluate the function
	assert(function);
  const Eigen::VectorXd& newOutput = function->Evaluate(newPt).at(0);

  // add the new point
  Add(newPt, newOutput);

  // return the result
  return newOutput;
}

unsigned int FlannCache::Add(Eigen::VectorXd const& input, Eigen::VectorXd const& result) {
  assert(input.size()==function->inputSizes(0));
  assert(result.size()==function->outputSizes(0));

	int cacheId = InCache(input);

	if(cacheId<0){
	  kdTree->add(input);
	  outputCache.push_back(result);

	  assert(kdTree->m_data.size()==outputCache.size());

		UpdateCentroid(input);

		return outputCache.size()-1;

	}else{
		return cacheId;
	}
}

void FlannCache::Remove(Eigen::VectorXd const& input) {
  // get the index of the point
  const int id = InCache(input);

  // the point is not in the cache ... nothing to remove
  if( id<0 ) { return; }

  // remove from output
  outputCache.erase(outputCache.begin()+id);

  // remove from input
  kdTree->m_data.erase(kdTree->m_data.begin()+id);
  kdTree->UpdateIndex();
}

size_t FlannCache::NearestNeighborIndex(Eigen::VectorXd const& point) const {
  // make sure we have enough
  assert(1<=Size());

  std::vector<size_t> indices;
  std::vector<double> squaredDists;
  std::tie(indices, squaredDists) = kdTree->query(point, 1);
  assert(indices.size()==1);
  assert(squaredDists.size()==1);

  return indices[0];
}

void FlannCache::NearestNeighbors(Eigen::VectorXd const& point,
                                  unsigned int const k,
                                  std::vector<Eigen::VectorXd>& neighbors,
                                  std::vector<Eigen::VectorXd>& result) const {
  // make sure we have enough
  assert(k<=Size());

  std::vector<size_t> indices;
  std::vector<double> squaredDists;
  std::tie(indices, squaredDists) = kdTree->query(point, k);
  assert(indices.size()==k);
  assert(squaredDists.size()==k);

  neighbors.resize(k);
  result.resize(k);
  for( unsigned int i=0; i<k; ++i ){
    neighbors.at(i) = kdTree->m_data.at(indices[i]);
    result.at(i) = outputCache.at(indices[i]);
  }
}

void FlannCache::NearestNeighbors(Eigen::VectorXd const& point,
                                  unsigned int const k,
                                  std::vector<Eigen::VectorXd>& neighbors) const {
  // make sure we have enough
  assert(k<=Size());

  std::vector<size_t> indices;
  std::vector<double> squaredDists;
  std::tie(indices, squaredDists) = kdTree->query(point, k);
  assert(indices.size()==k);
  assert(squaredDists.size()==k);

  neighbors.resize(k);
  for( unsigned int i=0; i<k; ++i ){ neighbors.at(i) = kdTree->m_data.at(indices[i]); }
}

unsigned int FlannCache::Size() const {
  // these two numbers should be the same unless we check the size after adding the input but before the model finishings running
  return std::min(kdTree->m_data.size(), outputCache.size());
}

std::vector<Eigen::VectorXd> FlannCache::Add(std::vector<Eigen::VectorXd> const& inputs) {
  std::vector<Eigen::VectorXd> results(inputs.size());

  for( unsigned int i=0; i<inputs.size(); ++i ) {
    // see if the point is already there
    const int index = InCache(inputs[i]);

    // add the point if is not already there
    results[i] = index<0? Add(inputs[i]) : outputCache.at(index);

    // make sure it got added
    assert(InCache(inputs[i])>=0);
  }

  assert(kdTree->m_data.size()==outputCache.size());

  return results;
}

void FlannCache::Add(std::vector<Eigen::VectorXd> const& inputs, std::vector<Eigen::VectorXd> const& results) {
  assert(inputs.size()==results.size());

  for( unsigned int i=0; i<inputs.size(); ++i ) {
    // add the point to cache (with result)
    Add(inputs[i], results[i]);

    // make sure it got added
    assert(InCache(inputs[i])>=0);
  }

  assert(kdTree->m_data.size()==outputCache.size());
}

const Eigen::VectorXd FlannCache::at(unsigned int const index) const {
  assert(index<kdTree->m_data.size());
  return kdTree->m_data[index];
}

Eigen::VectorXd FlannCache::at(unsigned int const index) {
  assert(index<kdTree->m_data.size());
  return kdTree->m_data[index];
}

Eigen::VectorXd const& FlannCache::OutputValue(unsigned int index) const{
	return outputCache.at(index);
}

void FlannCache::UpdateCentroid(Eigen::VectorXd const& point) {
	centroid = ((double)(Size()-1)*centroid+point)/(double)Size();
}

Eigen::VectorXd FlannCache::Centroid() const { return centroid; }

std::shared_ptr<ModPiece> FlannCache::Function() const { return function; }
