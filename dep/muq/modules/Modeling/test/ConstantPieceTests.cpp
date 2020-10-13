#include "gtest/gtest.h"

#include "MUQ/Modeling/ConstantPiece.h"

using namespace muq::Modeling;

/// A generic object for testing purposes
struct AnObject {
  AnObject(double const value) : value(value) {}

  const double value;
};

TEST(ConstantPiece, BasicTest) {
  // create an object
  auto obj = std::make_shared<AnObject>(2.0);

  // the outputs to the constant parameter
  std::vector<boost::any> outputs({1.0, (std::string)"string", obj});

  // create a constant parameter
  auto para = std::make_shared<ConstantPiece>(outputs);

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, 3);

  // evaluate the outputs
  auto outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 3);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs.at(0)), 1.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs.at(1)).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outs[2])->value, 2.0);
}

TEST(ConstantPiece, RecursiveConstructor) {
  // create an object
  auto obj = std::make_shared<AnObject>(2.0);

  // create a constant parameter
  auto para = std::make_shared<ConstantPiece>(1.0, (std::string)"string", obj);

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, 3);

  // evaluate the outputs
  auto outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 3);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs.at(0)), 1.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs.at(1)).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outs[2])->value, 2.0);
}

TEST(ConstantPiece, DefaultConstructor) {
  // create a constant parameter
  auto para = std::make_shared<ConstantPiece>();

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, -1);

  // evaluate, should return an empty vector
  auto outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 0);

  // create an object
  auto obj = std::make_shared<AnObject>(2.0);

  // the outputs to the constant parameter
  std::vector<boost::any> outputs({1.0, (std::string)"string", obj});

  // set new outputs
  para->SetOutputs(outputs);

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, -1);

  // evaluate the outputs
  outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 3);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs.at(0)), 1.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs.at(1)).compare("string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outs[2])->value, 2.0);

  // set new outputs
  para->SetOutputs(2.0, (std::string)"another string");

  // check the input/output number
  EXPECT_EQ(para->numInputs, 0);
  EXPECT_EQ(para->numOutputs, -1);

  // evaluate the outputs
  outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outs.at(0)), 2.0);
  EXPECT_TRUE(boost::any_cast<std::string>(outs.at(1)).compare("another string")==0);

  // set new outputs
  para->SetOutputs();

  // evaluate the outputs
  outs = para->Evaluate();

  // check that the outputs are what we expect
  EXPECT_EQ(outs.size(), 0);
}
