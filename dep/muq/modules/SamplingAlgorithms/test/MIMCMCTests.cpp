#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

#include <boost/property_tree/ptree.hpp>

#include <gtest/gtest.h>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;



class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> targetIn)
   : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,2), Eigen::VectorXi::Constant(1,2)),
     target(targetIn){}

  virtual ~MySamplingProblem() = default;


  virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> const& state, AbstractSamplingProblem::SampleType type) override {
    lastState = state;
    return target->Evaluate(state->state).at(0)(0);
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;

  std::shared_ptr<muq::Modeling::ModPiece> target;

};


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> Interpolate (std::shared_ptr<SamplingState> const& coarseProposal, std::shared_ptr<SamplingState> const& fineProposal) {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public MIComponentFactory {
public:
  virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
    0.6, 1.0;
    cov *= 20.0;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);
  }

  virtual std::shared_ptr<MultiIndex> FinestIndex() override {
    auto index = std::make_shared<MultiIndex>(2);
    index->SetValue(0, 2);
    index->SetValue(1, 2);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& index,
                                                        std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                           std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
    pt::ptree ptProposal;
    ptProposal.put("BlockIndex",0);
    int subsampling = 5;
    ptProposal.put("subsampling", subsampling);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
           0.6, 1.0;

    if (index->GetValue(0) == 0) {
      mu[0] *= 0.8;
    } else if (index->GetValue(0) == 1) {
      mu[0] *= 0.9;
    } else if (index->GetValue(0) == 2) {
      mu[0] *= 1.0;
    } else {
      std::cerr << "Sampling problem not defined!" << std::endl;
      assert (false);
    }
    if (index->GetValue(1) == 0) {
      mu[1] *= 0.8;
    } else if (index->GetValue(1) == 1) {
      mu[1] *= 0.9;
    } else if (index->GetValue(1) == 2) {
      mu[1] *= 1.0;
    } else {
      std::cerr << "Sampling problem not defined!" << std::endl;
      assert (false);
    }
    if (index->Max() == 0) {
      cov *= 2.0;
    } else if (index->Max() == 1) {
      cov *= 1.5;
    } else if (index->Max() == 2) {
      cov *= 1.0;
    } else {
      std::cerr << "Sampling problem not defined!" << std::endl;
      assert (false);
    }

    auto coarseTargetDensity = std::make_shared<Gaussian>(mu, cov)->AsDensity();
    return std::make_shared<MySamplingProblem>(coarseTargetDensity);
  }

  virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    return mu;
  }

};

TEST(MIMCMCTest, MIMCMC) {
  auto componentFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;
  pt.put("NumInitialSamples", 5e4); // number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); // estimator variance to be achieved by greedy algorithm

  MIMCMC mimcmc (pt, componentFactory);
  mimcmc.Run(Eigen::Vector2d(1.0, 2.0));
  mimcmc.Draw(false);

  auto mean = mimcmc.MeanQOI();

  EXPECT_NEAR(mean[0], 1.0, 0.25);
  EXPECT_NEAR(mean[1], 2.0, 0.25);

}

TEST(MIMCMCTest, SLMCMC)
{
  auto componentFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 5e3); // number of samples for single level

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.Run(Eigen::Vector2d(1.0, 2.0));

  auto mean = slmcmc.MeanQOI();

  EXPECT_NEAR(mean[0], 1.0, 0.1);
  EXPECT_NEAR(mean[1], 2.0, 0.1);

}
