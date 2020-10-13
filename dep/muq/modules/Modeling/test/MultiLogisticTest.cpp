#include "gtest/gtest.h"

#include "MUQ/Modeling/MultiLogisticLikelihood.h"
#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Modeling;
using namespace muq::Utilities;

TEST(MultiLogisticLikelihood, Evaluate)
{
  int numClasses = 3;

  int numData = 10;
  Eigen::VectorXi data = (numClasses*RandomGenerator::GetUniform(numData)).array().floor().cast<int>();

  auto piece = std::make_shared<MultiLogisticLikelihood>(numClasses,data);

  Eigen::VectorXd scores = RandomGenerator::GetNormal(numData*numClasses);
  Eigen::Map<Eigen::MatrixXd> scoreMat(scores.data(), numClasses, numData);

  double logLikely = piece->Evaluate(scores).at(0)(0);

  Eigen::MatrixXd logProbs = scoreMat.rowwise() - scoreMat.array().exp().colwise().sum().log().matrix();
  double trueLogLikely = 0.0;
  for(int i=0; i<numData;++i)
    trueLogLikely += logProbs(data(i),i);

  EXPECT_NEAR(trueLogLikely, logLikely, 1e-12);
}

TEST(MultiLogisticLikelihood, Derivatives)
{

  int numClasses = 2;

  int numData = 10;
  Eigen::VectorXi data = (numClasses*RandomGenerator::GetUniform(numData)).array().floor().cast<int>();

  auto piece = std::make_shared<MultiLogisticLikelihood>(numClasses,data);

  Eigen::VectorXd scores = RandomGenerator::GetNormal(numData*numClasses);

  Eigen::VectorXd sens = 0.5*Eigen::VectorXd::Ones(1);

  Eigen::VectorXd gradient = piece->Gradient(0,0, scores,sens);
  Eigen::VectorXd gradByFD = piece->GradientByFD(0,0,std::vector<Eigen::VectorXd>(1,scores), sens);

  for(int i=0; i<numData*numClasses; ++i)
    EXPECT_NEAR(gradByFD(i), gradient(i), 1e-3);

}
