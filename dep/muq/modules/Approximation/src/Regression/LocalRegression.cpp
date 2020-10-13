#include "MUQ/Approximation/Regression/LocalRegression.h"

#if MUQ_HAS_PARCER
#include <parcer/Eigen.h>
#endif

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Approximation;

LocalRegression::LocalRegression(std::shared_ptr<ModPiece> function, pt::ptree& pt) : ModPiece(function->inputSizes, function->outputSizes), kn(pt.get<unsigned int>("NumNeighbors")) {
  SetUp(function, pt);
}

#if MUQ_HAS_PARCER
LocalRegression::LocalRegression(std::shared_ptr<muq::Modeling::ModPiece> function, boost::property_tree::ptree& pt, std::shared_ptr<parcer::Communicator> comm) : ModPiece(function->inputSizes, function->outputSizes), kn(pt.get<unsigned int>("NumNeighbors")), comm(comm) /*tagSingle(comm->GetSize()+1), tagMulti(comm->GetSize()+2)*/ {
  SetUp(function, pt);
}

LocalRegression::~LocalRegression() {
  if( comm ) {
    // needs to be destroyed on all processors
    comm->Barrier();

    // clear any messages
    Probe();
  }
}
#endif

void LocalRegression::SetUp(std::shared_ptr<muq::Modeling::ModPiece> function, boost::property_tree::ptree& pt) {
  // can only have one input and output
  assert(inputSizes.size()==1);
  assert(outputSizes.size()==1);

  // create a cache of model evaluations
  cache = std::make_shared<FlannCache>(function);

  // create a regression object
  pt.put<std::string>("PolynomialBasis", pt.get<std::string>("PolynomialBasis", "Legendre")); // set default to Legendre
  pt.put<unsigned int>("Order", pt.get<unsigned int>("Order", 2)); // set default order to 2
  pt.put<double>("MaxPoisednessRadius", pt.get<double>("MaxPoisednessRadius", 1.0));
  pt.put<unsigned int>("InputSize", function->inputSizes(0));
  reg = std::make_shared<Regression>(pt);
}

void LocalRegression::ClearCache() {
  assert(cache);
  cache = std::make_shared<FlannCache>(cache->Function());
}

void LocalRegression::FitRegression(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  std::vector<Eigen::VectorXd> result;
  cache->NearestNeighbors(input, kn, neighbors, result);

  // fit the regression
  reg->Fit(neighbors, result, input);
}

void LocalRegression::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
#if MUQ_HAS_PARCER
  Probe();
#endif

  // fit the regressor
  FitRegression(inputs[0]);

  // evaluate the regressor
  outputs.resize(1);
  outputs[0] = (Eigen::VectorXd)boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(inputs[0].get()) [0]).col(0);
}

unsigned int LocalRegression::CacheSize() const {
  assert(cache);
  return cache->Size();
}

Eigen::VectorXd LocalRegression::CachePoint(unsigned int const index) const {
  assert(cache);
  return cache->at(index);
}

bool LocalRegression::InCache(Eigen::VectorXd const& point) const {
  assert(cache);
  return cache->InCache(point)>=0;
}

Eigen::VectorXd LocalRegression::Add(Eigen::VectorXd const& input) const {
  assert(cache);
  const Eigen::VectorXd result = cache->Add(input);

#if MUQ_HAS_PARCER
  if( comm ) {
    for( unsigned int i=0; i<comm->GetSize(); ++i ) {
      if( i==comm->GetRank() ) { continue; }

      comm->Send(std::pair<Eigen::VectorXd, Eigen::VectorXd>(input, result), i, tagSingle);
    }
    Probe();
  }
#endif

  return result;
}

void LocalRegression::Add(std::vector<Eigen::VectorXd> const& inputs) const {
  for( auto it : inputs ) { Add(it); }
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::PoisednessConstant(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  cache->NearestNeighbors(input, kn, neighbors);

  return PoisednessConstant(input, neighbors);
}

std::tuple<Eigen::VectorXd, double, unsigned int> LocalRegression::PoisednessConstant(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const {
  assert(neighbors.size()>=kn);
  assert(reg);
  std::pair<Eigen::VectorXd, double> lambda = reg->PoisednessConstant(neighbors, input, kn);

  double dist = 0.0;
  unsigned int index = 0;
  for( unsigned int i=0; i<kn; ++i ) {
    const double d = (input-neighbors[i]).norm();
    if( d<dist ) { dist = d; index = i; }
  }

  return std::tuple<Eigen::VectorXd, double, unsigned int>(lambda.first, lambda.second, index);
}

std::pair<double, double> LocalRegression::ErrorIndicator(Eigen::VectorXd const& input) const {
  // find the nearest neighbors
  std::vector<Eigen::VectorXd> neighbors;
  cache->NearestNeighbors(input, kn, neighbors);

  return ErrorIndicator(input, neighbors);
}

unsigned int LocalRegression::Order() const { return reg->order; }

std::pair<double, double> LocalRegression::ErrorIndicator(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const {
  assert(neighbors.size()==kn);

  // create a local factorial function (caution: may be problematic if n is sufficiently large)
  //std::function<int(int)> factorial = [&factorial](int const n) { return ((n==2 || n==1)? n : n*factorial(n-1)); };

  // compute the radius
  double delta = 0.0;
  for( auto n : neighbors) { delta = std::max(delta, (n-input).norm()); }
  /*Eigen::ArrayXd radius = Eigen::ArrayXd::Zero(input.size());
  for( auto n : neighbors) {
    //std::cout << n.transpose() << std::endl;
    radius = radius.max((n-input).array().abs());
  }*/

  return std::pair<double, double>(((double)reg->order+1.0)*std::log(delta), delta);
}

void LocalRegression::NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors) const {
  assert(cache);
  cache->NearestNeighbors(input, kn, neighbors);
}

void LocalRegression::NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& result) const {
  assert(cache);
  cache->NearestNeighbors(input, kn, neighbors, result);
}

Eigen::VectorXd LocalRegression::EvaluateRegressor(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors, std::vector<Eigen::VectorXd> const& result) const {
#if MUQ_HAS_PARCER
  Probe();
#endif

  // fit the regression
  reg->Fit(neighbors, result, input);

  // evaluate the regressor
  return (Eigen::VectorXd)boost::any_cast<Eigen::MatrixXd const&>(reg->Evaluate(input) [0]).col(0);
}

#if MUQ_HAS_PARCER
void LocalRegression::Probe() const {
  if( !comm ) { return; }
  //std::cout << "START PROBE!!! ";// /*<< comm->GetRank() */<< std::endl;

  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
  //  std::cout << "CHECKING " << i << " of " << comm->GetSize() << " rank: " << comm->GetRank();// << std::endl;
    //std::cout << "START RECIEVING " << i+1 << " of " << 2 << std::endl;
    if( i==comm->GetRank() ) { continue; }
    //std::cout << "START RECIEVING " << i << " of " << comm->GetSize() << " rank: " << comm->GetRank();// << std::endl;

    //{ // get single adds
      //parcer::RecvRequest recvReq;
      //bool enter = comm->Iprobe(i, tagSingle, recvReq);
      while( true ) {
        parcer::RecvRequest recvReq;
        if( !comm->Iprobe(i, tagSingle, recvReq) ) { break; }
        //std::cout << "PROBED ";
	      // get the point
        comm->Irecv(i, tagSingle, recvReq);
        //std::cout << "RECIEVED, proc: " << std::endl << std::flush;
        //std::cout << "RECIEVED, proc: " << comm->GetRank() << std::endl << std::flush;
        const std::pair<Eigen::VectorXd, Eigen::VectorXd> point = recvReq.GetObject<std::pair<Eigen::VectorXd, Eigen::VectorXd> >();
        //std::cout << std::flush;
        //auto point = recvReq.GetObject<std::pair<Eigen::VectorXd, Eigen::VectorXd> >();
        //std::cout << "GOT OBJECT!" << " rank: " << comm->GetRank() << " ";
        //std::cout << "GOT OBJECT!" << std::endl;

	      // add the point
	      cache->Add(point.first, point.second);

        //std::cout << "ADDED!" << " rank: " << comm->GetRank() << " ";

	      recvReq.Clear();

        //std::cout << "CLEARED!" << " rank: " << comm->GetRank() << " ";

        //enter = comm->Iprobe(i, tagSingle, recvReq);

      //  std::cout << "PROBED!" << " rank: " << comm->GetRank();
      }
    //}

    //std::cout << "DONE WITH THIS BIT" << " rank: " << comm->GetRank() << " ";
    //std::cout << "DONE RECIEVING " << i << " of " << comm->GetSize() << " ";// << std::endl;
  }

//  std::cout << "DONE WITH PROBE!!!" << std::endl;
}
#endif

Eigen::VectorXd LocalRegression::CacheCentroid() const {
  assert(cache);
  return cache->Centroid();
}
