#include <gtest/gtest.h>

#include "MUQ/Modeling/CwiseOperators/CwiseUnaryOperator.h"


using namespace muq::Modeling;


TEST(CwiseOperators, SimpleBase)
{
  unsigned int dim = 10;

  typedef CwiseUnaryOperator<std::exp, stan::math::exp> MyExpPiece;
  auto piece = std::make_shared<MyExpPiece>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::exp(input(i)), output(i));

  Eigen::MatrixXd jacobian = piece->Jacobian(0,0,input);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      if(i==j){
        EXPECT_DOUBLE_EQ(std::exp(input(i)), jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(0.0, jacobian(i,j));
      }
    }
  }

  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd gradient = piece->Gradient(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::exp(input(i)), gradient(i));

  Eigen::VectorXd jacAct = piece->ApplyJacobian(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::exp(input(i)), gradient(i));
}

TEST(CwiseOperators, ExpOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<ExpOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::exp(input(i)), output(i));

  Eigen::MatrixXd jacobian = piece->Jacobian(0,0,input);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      if(i==j){
        EXPECT_DOUBLE_EQ(std::exp(input(i)), jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(0.0, jacobian(i,j));
      }
    }
  }

  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd gradient = piece->Gradient(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::exp(input(i)), gradient(i));

  Eigen::VectorXd jacAct = piece->ApplyJacobian(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::exp(input(i)), gradient(i));
}

TEST(CwiseOperators, SinOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<SinOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::sin(input(i)), output(i));

  Eigen::MatrixXd jacobian = piece->Jacobian(0,0,input);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      if(i==j){
        EXPECT_DOUBLE_EQ(std::cos(input(i)), jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(0.0, jacobian(i,j));
      }
    }
  }

  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd gradient = piece->Gradient(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::cos(input(i)), gradient(i));

  Eigen::VectorXd jacAct = piece->ApplyJacobian(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i)*std::cos(input(i)), gradient(i));
}

TEST(CwiseOperators, CosOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<CosOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::cos(input(i)), output(i));

  Eigen::MatrixXd jacobian = piece->Jacobian(0,0,input);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      if(i==j){
        EXPECT_DOUBLE_EQ(-1.0*std::sin(input(i)), jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(0.0, jacobian(i,j));
      }
    }
  }

  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd gradient = piece->Gradient(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(-1.0*vec(i)*std::sin(input(i)), gradient(i));

  Eigen::VectorXd jacAct = piece->ApplyJacobian(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(-1.0*vec(i)*std::sin(input(i)), gradient(i));
}

TEST(CwiseOperators, AbsOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<AbsOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::abs(input(i)), output(i));

  Eigen::MatrixXd jacobian = piece->Jacobian(0,0,input);
  for(int j=0; j<dim; ++j){
    for(int i=0; i<dim; ++i){
      if(i==j){
        EXPECT_DOUBLE_EQ( (input(i)>0 ? 1.0 : -1.0), jacobian(i,j));
      }else{
        EXPECT_DOUBLE_EQ(0.0, jacobian(i,j));
      }
    }
  }

  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd gradient = piece->Gradient(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i) * (input(i)>0 ? 1.0 : -1.0), gradient(i));

  Eigen::VectorXd jacAct = piece->ApplyJacobian(0,0,input,vec);
  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(vec(i) * (input(i)>0 ? 1.0 : -1.0), gradient(i));
}

TEST(CwiseOperators, AcosOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<AcosOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::acos(input(i)), output(i));
}

TEST(CwiseOperators, CbrtOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<CbrtOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(std::cbrt(input(i)), output(i));
}

TEST(CwiseOperators, DigammaOperator)
{
  unsigned int dim = 10;

  auto piece = std::make_shared<DigammaOperator>(dim);

  Eigen::VectorXd input = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd output = piece->Evaluate(input).at(0);

  for(int i=0; i<dim; ++i)
    EXPECT_DOUBLE_EQ(stan::math::digamma(input(i)), output(i));
}
