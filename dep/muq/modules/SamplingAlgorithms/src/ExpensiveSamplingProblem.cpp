#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/SamplingAlgorithms/SamplingState.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, pt::ptree pt) : SamplingProblem(target), centroid(Eigen::VectorXd::Zero(target->inputSizes(0))) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));

  // must becaused after reg is created
  SetUp(pt);
}

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, Eigen::VectorXd const& centroid, pt::ptree pt) : SamplingProblem(target), centroid(centroid) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")));

  // must becaused after reg is created
  SetUp(pt);
}

#if MUQ_HAS_PARCER
ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm) : SamplingProblem(target), centroid(Eigen::VectorXd::Zero(target->inputSizes(0))) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")), comm);

  // must becaused after reg is created
  SetUp(pt);
}

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> const& target, Eigen::VectorXd const& centroid, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> comm) : SamplingProblem(target), centroid(centroid) {
  // create the local regressor
  reg = std::make_shared<LocalRegression>(target, pt.get_child(pt.get<std::string>("RegressionOptions")), comm);

  // must becaused after reg is created
  SetUp(pt);
}
#endif

void ExpensiveSamplingProblem::SetUp(boost::property_tree::ptree& pt) {
  assert(reg);

  // can only have one input
  assert(target->numInputs==1);

  beta = std::pair<double, double>(pt.get<double>("BetaScale", 0.0), -pt.get<double>("BetaExponent", RAND_MAX));
  assert(beta.second<0.0);

  tau0 = pt.get<double>("FirstLevelLength", 1.0);

  gamma = std::pair<double, double>(std::log(pt.get<double>("GammaScale", 1.0)), pt.get<double>("GammaExponent", 1.0));
  assert(gamma.second>0.0);

  lyapunovPara.first = pt.get<double>("LyapunovScale", 1.0);
  lyapunovPara.second = pt.get<double>("LyapunovExponent", 1.0);
  eta = pt.get<double>("TailCorrection", 0.0);
}

double ExpensiveSamplingProblem::LogDensity(unsigned int const step, std::shared_ptr<SamplingState> const& state, AbstractSamplingProblem::SampleType type) {
  std::vector<Eigen::VectorXd> neighbors, results;

  bool addThreshold = (type==AbstractSamplingProblem::SampleType::Proposed);

  double threshold = RefineSurrogate(step, state, neighbors, results);

  if( type==AbstractSamplingProblem::SampleType::Accepted ) {
    lastLyapunov = LogLyapunovFunction(state->state[0]);
    lastThreshold = threshold;
  }

  //double threshold = RefineSurrogate(step, state, neighbors, results);
  assert(neighbors.size()==results.size());
  assert(neighbors.size()>=reg->kn);

  // agument the threshold
  threshold *= (lastLyapunov<LogLyapunovFunction(state->state[0]) ? -1.0 : 1.0);

  // set cumulative refinement
  if( beta.first>0.0 ) { state->meta["cumulative random refinement"] = cumbeta; }
  state->meta["cumulative structural refinement"] = cumgamma;

  const double logtarg = reg->EvaluateRegressor(state->state[0],
                                std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn),
                                std::vector<Eigen::VectorXd>(results.begin(), results.begin()+reg->kn)) (0);

  return logtarg + (addThreshold ? eta*(lastThreshold+threshold)*state->state[0].norm()*std::fabs(logtarg) : 0.0);
}

double ExpensiveSamplingProblem::LogLyapunovFunction(Eigen::VectorXd const& x) const {
  /*(const Eigen::VectorXd diff = x-centroid;

  Eigen::VectorXd scale = Eigen::VectorXd::Constant(diff.size(), lyapunovScale);
  scale(0) = 1.0;
  scale(1) = 1.0;

  return diff.dot(scale.asDiagonal()*diff);*/

  //return lyapunovScale*(x-centroid).norm();
  return lyapunovPara.first*std::pow((x-centroid).norm(), lyapunovPara.second);
}

void ExpensiveSamplingProblem::CheckNumNeighbors(std::shared_ptr<SamplingState> const& state) {
  while( reg->CacheSize()<reg->kn ) {
    auto gauss = std::make_shared<muq::Modeling::Gaussian>(state->state[0], std::pow(std::exp(gamma.first), 1.0/(1.0+reg->Order()))*Eigen::VectorXd::Ones(state->state[0].size()));
    const Eigen::VectorXd& x = gauss->Sample();
    reg->Add(x);
    ++cumkappa;
  }
}

void ExpensiveSamplingProblem::CheckNeighbors(
  std::shared_ptr<SamplingState> const& state,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) const {
    reg->NearestNeighbors(state->state[0], neighbors, results);
    assert(neighbors.size()==results.size());
}

double ExpensiveSamplingProblem::ErrorThreshold(Eigen::VectorXd const& x) const {
  return LogLyapunovFunction(x) + gamma.first - gamma.second*std::log((double)level);
}

double ExpensiveSamplingProblem::RefineSurrogate(
  std::shared_ptr<SamplingState> const& state,
  double const radius,
  std::vector<Eigen::VectorXd>& neighbors,
  std::vector<Eigen::VectorXd>& results) {
    // compute the poisedness constant
    const std::tuple<Eigen::VectorXd, double, unsigned int>& lambda = reg->PoisednessConstant(state->state[0], neighbors);

    // if the point is already in the cache
    if( reg->InCache(std::get<0>(lambda)) ) {
      // choose a random point
      Eigen::VectorXd point = Eigen::VectorXd::Random(state->state[0].size());
      point *= RandomGenerator::GetUniform()*radius/point.norm();
      point += state->state[0];

      // find the closest point
      int index = 0;
      double dist = RAND_MAX;
      for( unsigned int i=0; i<reg->kn; ++i ) {
        double newdist = (point-state->state[0]).norm();
        if( newdist<dist ) { dist = newdist; index = i; }
      }

      assert(dist>0.0);
      assert(dist<=radius); // max is the diameter

      RefineSurrogate(point, index, neighbors, results);
    } else {
      RefineSurrogate(std::get<0>(lambda), std::get<2>(lambda), neighbors, results);
    }

  return std::get<1>(lambda);
}

double ExpensiveSamplingProblem::RefineSurrogate(unsigned int const step, std::shared_ptr<SamplingState> const& state, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) {
  // make sure we have enough points
  CheckNumNeighbors(state);

  // get the nearest nearest neighbors
  CheckNeighbors(state, neighbors, results);

  // get the error indicator (using the first kn neighbors)
  double error, radius;
  std::tie(error, radius) = reg->ErrorIndicator(state->state[0], neighbors);

  // random refinement
  if( RandomGenerator::GetUniform()<beta.first*std::pow((double)step, beta.second) ) {
    RefineSurrogate(state, radius, neighbors, results);
    ++cumbeta;

    // recompute the error indicator
    std::tie(error, radius) = reg->ErrorIndicator(state->state[0], std::vector<Eigen::VectorXd>(neighbors.begin(), neighbors.begin()+reg->kn));
  }

  // check to see if we should increment the level
  if( step>tau0*std::pow((double)level, 2.0*gamma.second) ) { ++level; }

  // compute (and store) the error threshold
  const double threshold = ErrorThreshold(state->state[0]);
  state->meta["error threshold"] = std::exp(threshold);
  state->meta["error indicator"] = std::exp(error);

  // structural refinement
  if( error>threshold ) {
    double lambda = RefineSurrogate(state, radius, neighbors, results);
    ++cumgamma;
  }

  return std::exp(threshold);
}

void ExpensiveSamplingProblem::RefineSurrogate(Eigen::VectorXd const& point, unsigned int const index, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) {
  assert(!reg->InCache(point));
  const Eigen::VectorXd& result = reg->Add(point);

  assert(neighbors.size()==results.size());
  assert(index<neighbors.size());

  neighbors.push_back(point);
  results.push_back(result);
  std::iter_swap(neighbors.end()-1, neighbors.begin()+index);
  std::iter_swap(results.end()-1, results.begin()+index);
}

unsigned int ExpensiveSamplingProblem::CacheSize() const { return reg->CacheSize(); }


void ExpensiveSamplingProblem::AddOptions(boost::property_tree::ptree & pt) const {
  pt.put("ReevaluateAcceptedDensity", true);
}
