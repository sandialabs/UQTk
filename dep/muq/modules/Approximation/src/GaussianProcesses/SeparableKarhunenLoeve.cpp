#include "MUQ/Approximation/GaussianProcesses/SeparableKarhunenLoeve.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

using namespace muq::Approximation;
using namespace muq::Utilities;

SeparableKarhunenLoeve::SeparableKarhunenLoeve(std::vector<std::shared_ptr<KernelBase>> kernels,
                                               std::vector<Eigen::MatrixXd> const& seedPts,
                                               std::vector<Eigen::VectorXd> const& seedWts,
                                               boost::property_tree::ptree options)
{

  const unsigned int numComps = kernels.size();

  // Check sizes of things
  assert(numComps == seedPts.size());
  assert(seedPts.size() == seedWts.size());

  const unsigned int inputDim = kernels.at(0).inputDim;

  for(int i=1; i<kernels.size(); ++i)
    assert(kernels.at(i).inputDim == inputDim);

  // Build a KL decomposition for each separable component
  components.resize(numComps);

  numModes = 1;

  Eigen::RowVectorXi allNumModes(numComps);

  for(int i=0; i<numComps; ++i){
    Eigen::MatrixXd newPts = Eigen::MatrixXd::Zero(inputDim, seedPts.at(i).cols());

    for(int j=0; j<kernels.at(i).dimInds.size(); ++j)
      newPts.row(kernels.at(i).dimInds.at(j)) = seedPts.at(i).row(j);

    components.at(i) = std::make_shared<KarhuneLoeveExpansion>(kernels.at(i), newPts, seedWts.at(i), options);

    allNumModes(i) = components.at(i)->NumModes();
    numModes *= components.at(i)->NumModes();
  }

  modeInds =  MultiIndexFactoryCreateFullTensor(allNumModes);

}

unsigned int SeparableKarhunenLoeve::NumModes()
{
  return numModes;
}

Eigen::MatrixXd SeparableKarhunenLoeve::GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts)
{
  std::vector<Eigen::MatrixXd> allModes(components.size());
  for(int i=0; i<components.size(); ++i)
    allModes.at(i) = components.at(i)->GetModes(pts);

  Eigen::MatrixXd output = Eigen::MatrixXd::Ones(pts.cols(), numModes);

  for(int i=0; i<modeInds->Size(); ++i){
    for(auto nzIt = modeInds->at(i)->GetNzBegin(); nzIt != modeInds->at(i)->GetNzEnd(); ++nzIt){
      output.col(i) *= allModes.at(nzIt->first).col(nzIt->second);
    }
  }

  return output;
}
