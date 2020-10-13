#include "MUQ/Modeling/SumPiece.h"

#include <gtest/gtest.h>

using namespace muq::Modeling;

TEST(Modeling, SumPiece)
{
  const int dim = 5;
  const int numInputs = 10;

  std::vector<Eigen::VectorXd> inputs(numInputs);
  Eigen::VectorXd trueSum = Eigen::VectorXd::Zero(dim);
  for(unsigned int i=0; i<numInputs; ++i){
    inputs.at(i) = Eigen::VectorXd::Random(dim);
    trueSum += inputs.at(i);
  }

  auto sum = std::make_shared<SumPiece>(dim, numInputs);

  Eigen::VectorXd output = sum->Evaluate(inputs).at(0);
  for(unsigned int d=0; d<dim; ++d)
    EXPECT_DOUBLE_EQ(trueSum(d), output(d));

  Eigen::MatrixXd jacobian = sum->Jacobian(0,1,inputs);
  for(unsigned int j=0; j<jacobian.cols(); ++j){
    for(unsigned int i=0; i<jacobian.rows(); ++i){
      if(i!=j){
        EXPECT_DOUBLE_EQ(0.0,jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(1.0,jacobian(i,j));
      }
    }
  }
}
