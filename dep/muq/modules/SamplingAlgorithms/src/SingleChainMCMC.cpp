#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"

#include <chrono>

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/StringUtilities.h"

#include "MUQ/SamplingAlgorithms/MarkovChain.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

SingleChainMCMC::SingleChainMCMC(pt::ptree pt,
                                 std::shared_ptr<AbstractSamplingProblem> const& problem) :
                 SamplingAlgorithm(std::make_shared<MarkovChain>()),
                 printLevel(pt.get("PrintLevel",3))
{
  Setup(pt, problem);
}

#if MUQ_HAS_PARCER
SingleChainMCMC::SingleChainMCMC(pt::ptree pt,
                                 std::shared_ptr<AbstractSamplingProblem> const& problem,
                                 std::shared_ptr<parcer::Communicator> const& comm) :
                 SamplingAlgorithm(std::make_shared<MarkovChain>(), comm),
                 printLevel(pt.get("PrintLevel",3))
{
  Setup(pt, problem);
}

SingleChainMCMC::SingleChainMCMC(pt::ptree pt,
                                 std::shared_ptr<parcer::Communicator> const& comm,
                                 std::vector<std::shared_ptr<TransitionKernel> > const& kernelsIn) :
                 SamplingAlgorithm(std::make_shared<MarkovChain>(), comm),
                 printLevel(pt.get("PrintLevel",3))
{
  Setup(pt, kernelsIn);
}

#endif

SingleChainMCMC::SingleChainMCMC(boost::property_tree::ptree pt,
                                 std::vector<std::shared_ptr<TransitionKernel> > const& kernelsIn) :
                SamplingAlgorithm(std::make_shared<MarkovChain>()),
                printLevel(pt.get("PrintLevel",3))
{
  Setup(pt,kernelsIn);
}

void SingleChainMCMC::Setup(pt::ptree pt,
                            std::vector<std::shared_ptr<TransitionKernel>> const& kernelsIn) {

  numSamps = pt.get<unsigned int>("NumSamples");
  burnIn = pt.get("BurnIn",0);

  kernels = kernelsIn;

  scheduler = std::make_shared<ThinScheduler>(pt);
  schedulerQOI = std::make_shared<ThinScheduler>(pt);
  assert(scheduler);
  assert(schedulerQOI);

#if MUQ_HAS_PARCER

  for(int i=0; i<kernels.size(); ++i)
    kernels.at(i)->SetCommunicator(comm);

#endif

}

void SingleChainMCMC::Setup(pt::ptree pt, std::shared_ptr<AbstractSamplingProblem> const& problem) {

  std::string kernelString = pt.get<std::string>("KernelList");

  std::vector<std::string> kernelNames = StringUtilities::Split(kernelString, ',');

  std::vector<std::shared_ptr<TransitionKernel>> kernelVec;
  unsigned int numBlocks = kernelNames.size();
  kernelVec.resize(numBlocks);

  // Add the block id to a child tree and construct a kernel for each block
  for(int i=0; i<numBlocks; ++i) {
    boost::property_tree::ptree subTree = pt.get_child(kernelNames.at(i));
    subTree.put("BlockIndex",i);

    problem->AddOptions(subTree);
    kernelVec.at(i) = TransitionKernel::Construct(subTree, problem);
  }

  Setup(pt, kernelVec);
}

void SingleChainMCMC::PrintStatus(std::string prefix, unsigned int currInd) const
{
  std::cout << prefix << int(std::floor(double((currInd - 1) * 100) / double(numSamps))) << "% Complete" << std::endl;
  if(printLevel>1){
    for(int blockInd=0; blockInd<kernels.size(); ++blockInd){
      std::cout << prefix << "  Block " << blockInd << ":\n";
      kernels.at(blockInd)->PrintStatus(prefix + "    ");
    }
  }
}

std::shared_ptr<SampleCollection> SingleChainMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
  if( !x0.empty() ) { SetState(x0); }

  // What is the next iteration that we want to print at
  const unsigned int printIncr = std::floor(numSamps / double(10));
  unsigned int nextPrintInd = printIncr;

  // Run until we've received the desired number of samples
  if(printLevel>0)
    std::cout << "Starting single chain MCMC sampler..." << std::endl;

  while(sampNum < numSamps)
  {
    // Should we print
    if(sampNum > nextPrintInd){
      if(printLevel>0){
        PrintStatus("  ", sampNum);
      }
      nextPrintInd += printIncr;
    }

    Sample();
  }


  if(printLevel>0){
    PrintStatus("  ", numSamps+1);
    std::cout << "Completed in " << totalTime << " seconds." << std::endl;
  }

  return samples;
}

void SingleChainMCMC::Sample() {
  auto startTime = std::chrono::high_resolution_clock::now();

  std::vector<std::shared_ptr<SamplingState> > newStates;

  // Loop through each parameter block
  for(int blockInd=0; blockInd<kernels.size(); ++blockInd){
    // kernel prestep
    kernels.at(blockInd)->PreStep(sampNum, prevState);

    // use the kernel to get the next state(s)
    newStates = kernels.at(blockInd)->Step(sampNum, prevState);

    // save when these samples where created
    const double now = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-startTime).count();
    for( auto it=newStates.begin(); it!=newStates.end(); ++it ) { (*it)->meta["time"] = now; }

    // kernel post-processing
    kernels.at(blockInd)->PostStep(sampNum, newStates);

    // add the new states to the SampleCollection (this also increments sampNum)
    prevState = SaveSamples(newStates, lastSavedState, sampNum);
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  totalTime += std::chrono::duration<double>(endTime - startTime).count();
}

std::shared_ptr<SamplingState> SingleChainMCMC::SaveSamples(std::vector<std::shared_ptr<SamplingState> > const& newStates, std::shared_ptr<SamplingState>& lastSavedState, unsigned int& sampNum) const {
  for( auto it : newStates ) {
    // save the sample, if we want to
    if( ShouldSave(sampNum) ) {
      samples->Add(it);
      if (it->HasMeta("QOI")) {
        std::shared_ptr<SamplingState> qoi = AnyCast(it->meta["QOI"]);
        QOIs->Add(qoi);
      }
    }

    // increment the number of samples and break of we hit the max. number
    if( ++sampNum>=numSamps ) { return it; }
  }

  return newStates.back();
}

bool SingleChainMCMC::ShouldSave(unsigned int const sampNum) const { return sampNum>=burnIn && scheduler->ShouldSave(sampNum); }

void SingleChainMCMC::SetState(std::vector<Eigen::VectorXd> const& x0) {
  prevState = std::make_shared<SamplingState>(x0);
  samples->Add(prevState);
}
