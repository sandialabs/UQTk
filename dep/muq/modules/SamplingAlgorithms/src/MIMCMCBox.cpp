#include "MUQ/SamplingAlgorithms/MIMCMCBox.h"

namespace muq {
  namespace SamplingAlgorithms {

    MIMCMCBox::MIMCMCBox(std::shared_ptr<MIComponentFactory> componentFactory, std::shared_ptr<MultiIndex> boxHighestIndex)
    : componentFactory(componentFactory),
    boxHighestIndex(boxHighestIndex)
    {
      pt::ptree ptChains;
      ptChains.put("NumSamples", 1e4); // number of MCMC steps
      pt::ptree ptBlockID;
      ptBlockID.put("BlockIndex",0);

      const auto rootIndex = std::make_shared<MultiIndex>(boxHighestIndex->GetLength());

      // Set up root index sampling
      auto coarse_problem = componentFactory->SamplingProblem(rootIndex);
      auto proposal_coarse = componentFactory->Proposal(rootIndex, coarse_problem);

      std::vector<std::shared_ptr<TransitionKernel>> coarse_kernels(1);
      coarse_kernels[0] = std::make_shared<MHKernel>(ptBlockID,coarse_problem,proposal_coarse);

      Eigen::VectorXd startPtCoarse = componentFactory->StartingPoint(rootIndex);
      auto coarse_chain = std::make_shared<SingleChainMCMC>(ptChains,coarse_kernels);
      coarse_chain->SetState(startPtCoarse);

      // Construct path to lowest index of box
      boxLowestIndex = MultiIndex::Copy(boxHighestIndex);
      --(*boxLowestIndex);
      std::shared_ptr<MultiIndexSet> rootPath = CreateRootPath(boxLowestIndex);
      for (int i = rootPath->Size()-2; i >= 0; i--) {

        tailChains.push_back(coarse_chain);

        std::shared_ptr<MultiIndex> index = (*rootPath)[i];

        auto problem = componentFactory->SamplingProblem(index);
        auto proposal = componentFactory->Proposal(index, problem);
        auto coarse_proposal = componentFactory->CoarseProposal(index, coarse_problem, coarse_chain);
        auto proposalInterpolation = componentFactory->Interpolation(index);
        auto startingPoint = componentFactory->StartingPoint(index);

        std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
        kernels[0] = std::make_shared<MIKernel>(ptBlockID,problem,coarse_problem,proposal,coarse_proposal,proposalInterpolation,coarse_chain);

        auto chain = std::make_shared<SingleChainMCMC>(ptChains,kernels);
        chain->SetState(startingPoint);

        coarse_problem = problem;
        coarse_chain = chain;
      }

      finestProblem = coarse_problem;

      std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);

      // Set up Multiindex box
      boxIndices = MultiIndexFactory::CreateFullTensor(boxSize->GetVector());
      boxChains.resize(boxIndices->Size());

      for (int i = 0; i < boxIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];

        if (boxIndex->Max() == 0) {
          boxChains[boxIndices->MultiToIndex(boxIndex)] = coarse_chain;
          continue;
        }

        std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);

        auto problem = componentFactory->SamplingProblem(index);
        auto proposal = componentFactory->Proposal(index, problem);
        auto coarse_proposal = componentFactory->CoarseProposal(index, coarse_problem, coarse_chain);
        auto proposalInterpolation = componentFactory->Interpolation(index);
        auto startingPoint = componentFactory->StartingPoint(index);

        std::vector<std::shared_ptr<TransitionKernel>> kernels(1);
        kernels[0] = std::make_shared<MIKernel>(ptBlockID,problem,coarse_problem,proposal,coarse_proposal,proposalInterpolation,coarse_chain);

        auto chain = std::make_shared<SingleChainMCMC>(ptChains,kernels);
        chain->SetState(startingPoint);

        boxChains[boxIndices->MultiToIndex(boxIndex)] = chain;

        if (boxIndex->Max() == 1)
          finestProblem = problem;
      }
    }

    std::shared_ptr<AbstractSamplingProblem> MIMCMCBox::GetFinestProblem() {
      return finestProblem;
    }


    void MIMCMCBox::Sample() {
      for (int i = 0; i < boxIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
        auto chain = boxChains[boxIndices->MultiToIndex(boxIndex)];
        chain->Sample();
      }
    }

    Eigen::VectorXd MIMCMCBox::MeanQOI() {
      Eigen::VectorXd sampMean = Eigen::VectorXd::Zero(GetFinestProblem()->blockSizesQOI.sum());

      for (int i = 0; i < boxIndices->Size(); i++) {
        std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
        auto chain = boxChains[boxIndices->MultiToIndex(boxIndex)];
        auto samps = chain->GetQOIs();

        std::shared_ptr<MultiIndex> index = std::make_shared<MultiIndex>(*boxLowestIndex + *boxIndex);
        auto indexDiffFromTop = std::make_shared<MultiIndex>(*boxHighestIndex - *index);

        if (indexDiffFromTop->Sum() % 2 == 0) {
          sampMean += samps->Mean();
        } else {
          sampMean -= samps->Mean();
        }
      }
      return sampMean;
    }

    void MIMCMCBox::DrawChain(std::shared_ptr<SingleChainMCMC> chain, std::string chainid, std::ofstream& graphfile) const {
      graphfile << "subgraph cluster_" << chainid << " {" << std::endl;
      graphfile << "label=\"Chain " << chainid << "\"" << std::endl;
      for (int s = 0; s < chain->GetSamples()->size(); s++) {
        std::shared_ptr<SamplingState> sample = chain->GetSamples()->at(s);
        std::string nodeid = "\"s" + chainid + "node" + std::to_string(s) + "\"";
        sample->meta["gvizid"] = nodeid;

        double logTarget = AnyCast(sample->meta["LogTarget"]);
        graphfile << nodeid << " [label=\""
        << s << " - " << sample->weight
        << " (L=" << logTarget << ")";

        if (sample->HasMeta("QOI"))
          graphfile << " QOI";

        graphfile << "\"]" << std::endl;
      }
      graphfile << "}" << std::endl;
      for (int s = 0; s < chain->GetSamples()->size(); s++) {
        std::shared_ptr<SamplingState> sample = chain->GetSamples()->at(s);
        std::string nodeid = "\"s" + chainid + "node" + std::to_string(s) + "\"";

        if (s < chain->GetSamples()->size() - 1)
          graphfile << nodeid << " -> " << "\"s" << chainid << "node" << s+1 << "\"" << std::endl;

        if (sample->HasMeta("coarseSample")) {
          std::shared_ptr<SamplingState> coarseSample = AnyCast(sample->meta["coarseSample"]);
          if (coarseSample->HasMeta("gvizid")) {
            std::string coarseid = AnyCast(coarseSample->meta["gvizid"]);
            graphfile << nodeid << " -> " << coarseid << std::endl;
          } else {
            std::cout << "no gvizid!" << std::endl;
          }
        }
      }
    }

    void MIMCMCBox::Draw(std::ofstream& graphfile, bool drawSamples) const {

      if (drawSamples) {
        for (int i = 0; i < tailChains.size(); i++) {
          std::string chainid = "box" + boxHighestIndex->ToString() + "_tail" + std::to_string(i) + "";

          DrawChain (tailChains[i], chainid, graphfile);
        }

        for (int i = 0; i < boxIndices->Size(); i++) {
          std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
          std::shared_ptr<SingleChainMCMC> singleChain = boxChains[boxIndices->MultiToIndex(boxIndex)];

          std::string chainid = "box" + boxHighestIndex->ToString() + "_node" + boxIndex->ToString() + "";
          DrawChain (singleChain, chainid, graphfile);
        }
      } else {
        std::string previd = "";
        for (int i = 0; i < tailChains.size(); i++) {
          std::string chainid = "\"box" + boxHighestIndex->ToString() + "_tail" + std::to_string(i) + "\"";
          if (i > 0)
            graphfile << chainid << " -> " << previd << std::endl;
          else
            graphfile << chainid << std::endl;
          previd = chainid;
        }

        for (int i = 0; i < boxIndices->Size(); i++) {
          std::shared_ptr<MultiIndex> boxIndex = (*boxIndices)[i];
          std::shared_ptr<SingleChainMCMC> singleChain = boxChains[boxIndices->MultiToIndex(boxIndex)];

          std::string chainid = "\"box" + boxHighestIndex->ToString() + "_node" + boxIndex->ToString() + "\"";

          if (previd != "")
            graphfile << chainid << " -> " << previd << std::endl;
          else
            graphfile << chainid << std::endl;
          if (i == 0)
            previd = chainid;
        }
      }
    }

    std::shared_ptr<SingleChainMCMC> MIMCMCBox::FinestChain() {
      std::shared_ptr<MultiIndex> boxSize = std::make_shared<MultiIndex>(*boxHighestIndex - *boxLowestIndex);
      return boxChains[boxIndices->MultiToIndex(boxSize)];
    }

    std::shared_ptr<MultiIndexSet> MIMCMCBox::CreateRootPath(std::shared_ptr<MultiIndex> index) {

      // create an empy multiindex set
      std::shared_ptr<MultiIndexLimiter> limiter = std::make_shared<NoLimiter>();
      std::shared_ptr<MultiIndexSet> output = std::make_shared<MultiIndexSet>(index->GetLength(),limiter);

      // Always go down one step in direction of largest index value until reaching root index
      std::shared_ptr<MultiIndex> currentIndex = index;
      output->AddActive(currentIndex);

      while (true) {

        int maxCoeffId;
        int maxEntry = currentIndex->GetVector().maxCoeff(&maxCoeffId);
        if (maxEntry == 0)
          break;

        std::shared_ptr<MultiIndex> nextIndex = MultiIndex::Copy(currentIndex);
        nextIndex->SetValue(maxCoeffId, maxEntry - 1);
        output->AddActive(nextIndex);

        currentIndex = nextIndex;
      }

      return output;
    }

  }
}
