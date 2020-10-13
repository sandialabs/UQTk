#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

using namespace muq::SamplingAlgorithms;

SamplingAlgorithm::SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samplesIn, std::shared_ptr<SampleCollection> const& QOIsIn) : samples(samplesIn), QOIs(QOIsIn) {}

SamplingAlgorithm::SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samplesIn) : SamplingAlgorithm(samplesIn, std::make_shared<SampleCollection>()) {}

#if MUQ_HAS_PARCER
SamplingAlgorithm::SamplingAlgorithm(std::shared_ptr<SampleCollection> const& samplesIn, std::shared_ptr<parcer::Communicator> const& comm) : samples(samplesIn), comm(comm) {}
#endif

std::shared_ptr<SampleCollection> SamplingAlgorithm::GetSamples() const { return samples; }

std::shared_ptr<SampleCollection> SamplingAlgorithm::GetQOIs() const { return QOIs; }

void SamplingAlgorithm::SetState(std::vector<Eigen::VectorXd> const& x0) {}

std::shared_ptr<SampleCollection> SamplingAlgorithm::Run(std::vector<Eigen::VectorXd> const& x0) { return RunImpl(x0); };

#if MUQ_HAS_PARCER
std::shared_ptr<parcer::Communicator> SamplingAlgorithm::GetCommunicator() const { return comm; }
#endif
