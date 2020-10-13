#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/WorkPiece.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;

REGISTER_MCMC_PROPOSAL(CrankNicolsonProposal)

CrankNicolsonProposal::CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                                             std::shared_ptr<AbstractSamplingProblem> prob,
                                             std::shared_ptr<GaussianBase>                priorIn) : MCMCProposal(pt,prob),
                                                                                                     beta(pt.get("Beta",0.5)),
                                                                                                     priorDist(priorIn)
{
}

CrankNicolsonProposal::CrankNicolsonProposal(boost::property_tree::ptree       const& pt,
                                             std::shared_ptr<AbstractSamplingProblem> prob) : MCMCProposal(pt,prob),
                                                                                              beta(pt.get("Beta",0.5))
{
  ExtractPrior(prob, pt.get<std::string>("PriorNode"));
}


std::shared_ptr<SamplingState> CrankNicolsonProposal::Sample(std::shared_ptr<SamplingState> currentState)
{
  // the mean of the proposal is the current point
  std::vector<Eigen::VectorXd> props = currentState->state;

  std::vector<Eigen::VectorXd> hypers = GetPriorInputs(currentState->state);

  Eigen::VectorXd priorSamp = priorDist->Sample(hypers);

  props.at(blockInd) = priorDist->GetMean() + sqrt(1.0-beta*beta)*(currentState->state.at(blockInd)-priorDist->GetMean()) + beta*(priorSamp - priorDist->GetMean());

  // store the new state in the output
  return std::make_shared<SamplingState>(props, 1.0);
}

double CrankNicolsonProposal::LogDensity(std::shared_ptr<SamplingState> currState,
                                         std::shared_ptr<SamplingState> propState)
{
  std::vector<Eigen::VectorXd> hypers = GetPriorInputs(currState->state);
  if(hypers.size()>0)
    priorDist->ResetHyperparameters(WorkPiece::ToRefVector(hypers));

  Eigen::VectorXd diff = (propState->state.at(blockInd) - priorDist->GetMean() - sqrt(1.0-beta*beta)*(currState->state.at(blockInd)-priorDist->GetMean()))/beta;

  hypers.insert(hypers.begin(), (diff + priorDist->GetMean()).eval());
  return priorDist->LogDensity(hypers);
}

std::vector<Eigen::VectorXd> CrankNicolsonProposal::GetPriorInputs(std::vector<Eigen::VectorXd> const& currState)
{
  std::vector<Eigen::VectorXd> hyperParams;

  if(priorMeanModel){
    ref_vector<Eigen::VectorXd> meanIns;
    for(int i=0; i<priorMeanInds.size(); ++i)
      meanIns.push_back( std::cref( currState.at(priorMeanInds.at(i)) ) );

    hyperParams.push_back(priorMeanModel->Evaluate(meanIns).at(0));
  }

  if(priorCovModel){
    ref_vector<Eigen::VectorXd> covIns;
    for(int i=0; i<priorCovInds.size(); ++i)
      covIns.push_back( std::cref( currState.at(priorCovInds.at(i)) ) );

    hyperParams.push_back(priorCovModel->Evaluate(covIns).at(0));
  }

  return hyperParams;
}


void CrankNicolsonProposal::ExtractPrior(std::shared_ptr<AbstractSamplingProblem> prob,
                                         std::string                              nodeName)
{
  // Cast the abstract base class into a sampling problem
  std::shared_ptr<SamplingProblem> prob2 = std::dynamic_pointer_cast<SamplingProblem>(prob);
  assert(prob2);

  // From the sampling problem, extract the ModPiece and try to cast it to a ModGraphPiece
  std::shared_ptr<ModPiece> targetDens = prob2->GetDistribution();
  std::shared_ptr<ModGraphPiece> targetDens2 = std::dynamic_pointer_cast<ModGraphPiece>(targetDens);
  assert(targetDens2);

  // Get the graph
  auto graph = targetDens2->GetGraph();

  // Get the prior piece corresponding to the Gaussian name
  auto priorPiece = graph->GetPiece(nodeName);
  assert(priorPiece);

  // Get the prior distribution
  std::shared_ptr<Density> priorDens = std::dynamic_pointer_cast<Density>(priorPiece);
  assert(priorDens);

  priorDist = std::dynamic_pointer_cast<GaussianBase>(priorDens->GetDistribution());
  assert(priorDist);

  // Check to see if the prior has a mean or covariance input.
  auto gaussPrior = std::dynamic_pointer_cast<Gaussian>(priorDist);
  if(gaussPrior==nullptr)
    return;

  Gaussian::InputMask inputTypes = gaussPrior->GetInputTypes();

  if(inputTypes == Gaussian::None)
    return;

  auto newGraph = graph->Clone();
  auto constPieces = targetDens2->GetConstantPieces();
  for(int i=0; i<constPieces.size(); ++i)
    newGraph->RemoveNode( newGraph->GetName(constPieces.at(i)) );

  // Check if there is a mean input
  if(inputTypes & Gaussian::Mean){
    std::string meanInput = graph->GetParent(nodeName, 1);
    if(newGraph->HasNode(meanInput)){
      priorMeanModel = newGraph->CreateModPiece(meanInput);
      priorMeanInds = targetDens2->MatchInputs(std::dynamic_pointer_cast<ModGraphPiece>(priorMeanModel));
    }else{
      priorMeanModel = std::make_shared<IdentityOperator>(priorDens->inputSizes(1));
      auto iter = std::find(constPieces.begin(), constPieces.end(), graph->GetPiece(meanInput));
      priorMeanInds.push_back( std::distance(constPieces.begin(),iter));
    }
  }

  // Check if there is a covariance or precision input
  if((inputTypes & Gaussian::DiagCovariance) || (inputTypes & Gaussian::FullCovariance) || (inputTypes & Gaussian::DiagPrecision)||(inputTypes & Gaussian::FullPrecision)){
    int covInd = (inputTypes & Gaussian::Mean) ? 2 : 1;
    priorUsesCov = (inputTypes & Gaussian::DiagCovariance) || (inputTypes & Gaussian::FullCovariance);

    std::string covInput = graph->GetParent(nodeName, covInd);
    if(newGraph->HasNode(covInput)){
      priorCovModel = newGraph->CreateModPiece(covInput);
      priorCovInds = targetDens2->MatchInputs(std::dynamic_pointer_cast<ModGraphPiece>(priorCovModel));
    }else{
      priorCovModel = std::make_shared<IdentityOperator>(priorDens->inputSizes(covInd));
      auto iter = std::find(constPieces.begin(), constPieces.end(), graph->GetPiece(covInput));
      priorCovInds.push_back( std::distance(constPieces.begin(),iter));
    }
  }

  // if(priorDens->inputSizes.size()>1){
  //
  //   // Create a new graph that does not have the constant pieces
  //   auto newGraph = graph->Clone();
  //   auto constPieces = targetDens2->GetConstantPieces();
  //   for(int i=0; i<constPieces.size(); ++i)
  //     newGraph->RemoveNode( newGraph->GetName(constPieces.at(i)) );
  //
  //   std::vector<std::string> gaussInputs  = graph->GetParents(nodeName);
  //   inputPieces.resize(gaussInputs.size()-1);
  //   blockInds.resize(gaussInputs.size()-1);
  //
  //   std::vector<std::vector<int>> blockInds(gaussInputs.size()-1);
  //
  //   for(int i=0; i<inputPieces.size(); ++i){
  //
  //     if(newGraph->HasNode(gaussInputs.at(i+1))){
  //       inputPieces.at(i) = newGraph->CreateModPiece(gaussInputs.at(i+1));
  //       blockInds.at(i) = targetDens2->MatchInputs(std::dynamic_pointer_cast<ModGraphPiece>(inputPieces.at(i)));
  //     }else{
  //       inputPieces.at(i) = std::make_shared<IdentityOperator>(priorDens->inputSizes(i));
  //       auto iter = std::find(constPieces.begin(), constPieces.end(), graph->GetPiece(gaussInputs.at(i+1)));
  //       blockInds.at(i).push_back( std::distance(constPieces.begin(),iter));
  //     }
  //   }
  //
  // }

}
