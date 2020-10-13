#include "MUQ/SamplingAlgorithms/ISKernel.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

REGISTER_TRANSITION_KERNEL(ISKernel)

ISKernel::ISKernel(pt::ptree const& pt, std::shared_ptr<SamplingProblem> problem) : TransitionKernel(pt, problem) {}

ISKernel::~ISKernel() {}

void ISKernel::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  const boost::any state = problem->SampleBiasingDistribution(inputs);

  ref_vector<boost::any> dist_inputs(inputs.begin(), inputs.end());
  dist_inputs.insert(dist_inputs.begin(), std::cref(state));

  const double logTarget = problem->EvaluateLogTarget(dist_inputs);
  const double logBias = problem->EvaluateLogBiasingDistribution(dist_inputs);

  outputs.resize(1);
  outputs[0] = std::make_shared<SamplingState>(state, std::exp(logTarget-logBias));
}
