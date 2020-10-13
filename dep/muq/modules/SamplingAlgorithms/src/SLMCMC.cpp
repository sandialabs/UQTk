#include "MUQ/SamplingAlgorithms/SLMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    SLMCMC::SLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
     : SamplingAlgorithm(std::shared_ptr<SampleCollection>(), std::shared_ptr<SampleCollection>()),
       componentFactory(componentFactory)
    {
      auto index = componentFactory->FinestIndex();

      auto problem = componentFactory->SamplingProblem(index);
      auto proposal = componentFactory->Proposal(index, problem);

      pt::ptree ptBlockID;
      ptBlockID.put("BlockIndex",0);
      std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
      kernels[0] = std::make_shared<MHKernel>(ptBlockID,problem,proposal);

      Eigen::VectorXd startingPoint = componentFactory->StartingPoint(index);

      coarse_chain = std::make_shared<SingleChainMCMC>(pt,kernels);
      coarse_chain->SetState(startingPoint);
    }

    std::shared_ptr<SampleCollection> SLMCMC::GetSamples() const {
      return nullptr;
    }
    std::shared_ptr<SampleCollection> SLMCMC::GetQOIs() const {
      return nullptr;
    }

    std::shared_ptr<SampleCollection> SLMCMC::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
      coarse_chain->Run();
      return nullptr;
    }

    Eigen::VectorXd SLMCMC::MeanQOI() {
      return coarse_chain->GetQOIs()->Mean();
    }

    Eigen::VectorXd SLMCMC::MeanParameter() {
      auto samps = coarse_chain->GetSamples();
      return samps->Mean();
    }

  }
}
