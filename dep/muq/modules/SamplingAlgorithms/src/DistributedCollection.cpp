#include "MUQ/SamplingAlgorithms/DistributedCollection.h"

#if MUQ_HAS_MPI

#if MUQ_HAS_PARCER
#include "parcer/Eigen.h"
#endif

using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

DistributedCollection::DistributedCollection(std::shared_ptr<SampleCollection> collection, std::shared_ptr<parcer::Communicator> comm) : SampleCollection(), collection(collection), comm(comm) {}

void DistributedCollection::Add(std::shared_ptr<SamplingState> newSamp) {
  collection->Add(newSamp);
}

std::shared_ptr<SamplingState> DistributedCollection::at(unsigned i) { return GlobalAt(i); }

const std::shared_ptr<SamplingState> DistributedCollection::at(unsigned i) const { return GlobalAt(i); }

std::shared_ptr<SamplingState> DistributedCollection::LocalAt(unsigned i) {
  assert(i<collection->size());
  return collection->at(i);
}

const std::shared_ptr<SamplingState> DistributedCollection::LocalAt(unsigned i) const {
  assert(i<collection->size());
  return collection->at(i);
}

std::shared_ptr<SamplingState> DistributedCollection::GlobalAt(unsigned i) {
  assert(i<GlobalSize());

  std::shared_ptr<SamplingState> state = nullptr;

  int size = 0;
  for( unsigned int j=0; j<comm->GetSize(); ++j ) {
    int localSize = j==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, j);

    if( i<size+localSize ) {
      state = j==comm->GetRank() ? LocalAt(i-size) : nullptr;
      comm->Bcast(state, j);
      break;
    }

    size += localSize;
  }

  return state;
}

const std::shared_ptr<SamplingState> DistributedCollection::GlobalAt(unsigned i) const {
  assert(i<GlobalSize());

  std::shared_ptr<SamplingState> state = nullptr;

  int size = 0;
  for( unsigned int j=0; j<comm->GetSize(); ++j ) {
    int localSize = j==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, j);

    if( i<size+localSize ) {
      state = j==comm->GetRank() ? LocalAt(i-size) : nullptr;
      comm->Bcast(state, j);
      break;
    }

    size += localSize;
  }

  return state;
}

unsigned int DistributedCollection::LocalSize() const { return collection->size(); }

unsigned int DistributedCollection::GlobalSize() const {
  int size = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    size += localSize;
  }

  return size;
}

unsigned int DistributedCollection::size() const { return GlobalSize(); }

Eigen::VectorXd DistributedCollection::LocalCentralMoment(unsigned order, int blockDim) const { return collection->CentralMoment(order, GlobalMean(blockDim), blockDim); }

Eigen::VectorXd DistributedCollection::GlobalCentralMoment(unsigned order, int blockDim) const { return GlobalEigenMean(LocalCentralMoment(order, blockDim)); }

Eigen::VectorXd DistributedCollection::CentralMoment(unsigned order, int blockDim) const { return GlobalCentralMoment(order, blockDim); }

Eigen::VectorXd DistributedCollection::LocalMean(int blockDim) const { return collection->Mean(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalMean(int blockDim) const { return GlobalEigenMean(LocalMean(blockDim)); }

Eigen::VectorXd DistributedCollection::Mean(int blockDim) const { return GlobalMean(blockDim); }

Eigen::VectorXd DistributedCollection::LocalVariance(int blockDim) const { return collection->Variance(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalVariance(int blockDim) const { return GlobalEigenMean(LocalVariance(blockDim)); }

Eigen::VectorXd DistributedCollection::Variance(int blockDim) const { return GlobalVariance(blockDim); }

Eigen::MatrixXd DistributedCollection::LocalCovariance(int blockDim) const { return collection->Covariance(GlobalMean(blockDim), blockDim); }

Eigen::MatrixXd DistributedCollection::GlobalCovariance(int blockDim) const { return GlobalEigenMean(LocalCovariance(blockDim)); }

Eigen::MatrixXd DistributedCollection::Covariance(int blockDim) const { return GlobalCovariance(blockDim); }

Eigen::VectorXd DistributedCollection::LocalESS(int blockDim) const { return collection->ESS(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalESS(int blockDim) const {
  const Eigen::VectorXd& local = LocalESS(blockDim);
  Eigen::VectorXd global = Eigen::VectorXd::Zero(local.size());

  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    Eigen::VectorXd l(local.size());
    if( comm->GetRank()==i ) { l = local; }
    comm->Bcast(l, i);

    global += l;
  }

  return global;
}

Eigen::VectorXd DistributedCollection::ESS(int blockDim) const { return GlobalESS(blockDim); }

Eigen::MatrixXd DistributedCollection::AsLocalMatrix(int blockDim) const { return collection->AsMatrix(blockDim); }

Eigen::MatrixXd DistributedCollection::AsGlobalMatrix(int blockDim) const {
  const Eigen::MatrixXd& local = AsLocalMatrix(blockDim);
  Eigen::MatrixXd global(local.rows(), GlobalSize());

  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    Eigen::MatrixXd l(local.rows(), localSize);
    if( comm->GetRank()==i ) { l = local; }

    comm->Bcast(l, i);
    global.block(0, numSamps, local.rows(), localSize) = l;

    numSamps += localSize;
  }

  return global;
}

Eigen::MatrixXd DistributedCollection::AsMatrix(int blockDim) const { return AsGlobalMatrix(blockDim); }

Eigen::VectorXd DistributedCollection::LocalWeights() const { return collection->Weights(); }

Eigen::VectorXd DistributedCollection::GlobalWeights() const {
  const Eigen::VectorXd& local = LocalWeights();
  Eigen::VectorXd global = Eigen::VectorXd::Constant(GlobalSize(), std::numeric_limits<double>::quiet_NaN());

  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    Eigen::VectorXd l(localSize);
    if( comm->GetRank()==i ) { l = local; }

    comm->Bcast(l, i);
    global.segment(numSamps, localSize) = l;

    numSamps += localSize;
  }

  return global;
}

Eigen::VectorXd DistributedCollection::Weights() const {
  return GlobalWeights();
}

void DistributedCollection::WriteToFile(std::string const& filename, std::string const& dataset) const {
  collection->WriteToFile(filename, dataset);
}

Eigen::VectorXd DistributedCollection::GlobalExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains) const {
  const Eigen::VectorXd& local = LocalExpectedValue(f, metains);
  Eigen::VectorXd global = Eigen::VectorXd::Zero(f->outputSizes(0));

  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    Eigen::VectorXd l(f->outputSizes(0));
    if( comm->GetRank()==i ) { l = local; }

    comm->Bcast(l, i);
    global += l;
  }

  return global/(double)comm->GetSize();
}

Eigen::VectorXd DistributedCollection::LocalExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains) const {
  return collection->ExpectedValue(f, metains);
}

Eigen::VectorXd DistributedCollection::ExpectedValue(std::shared_ptr<muq::Modeling::ModPiece> const& f, std::vector<std::string> const& metains) const {
  return GlobalExpectedValue(f, metains);
}

#endif // end MUQ_HAS_MPI
