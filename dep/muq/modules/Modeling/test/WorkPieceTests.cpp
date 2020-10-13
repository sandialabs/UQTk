#include "gtest/gtest.h"

#include "WorkPieceTestClasses.h"

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class WorkPieceTests : public::testing::Test {
public:

  /// Default constructor
  inline WorkPieceTests() {}

  /// Default destructor
  inline virtual ~WorkPieceTests() {}

  /// A string type object to input
  const std::string s = "a string";

  /// A double type object to input
  const double a = 2.0;

  /// An unsigned int type object to input
  const unsigned int b = 3;

private:
};

TEST_F(WorkPieceTests, Unfixed) {
  // create the test WorkPiece
  auto test = std::make_shared<UnfixedMod>();

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // evaluate with zero inputs and three outputs
  std::vector<boost::any> outputs = test->Evaluate();

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 3);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare("hello!")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 3.0);
  EXPECT_EQ(boost::any_cast<int>(outputs[2]), 6);

  // evaluate with one input and one output
  outputs = test->Evaluate(s);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);

  // evaluate with two inputs and zero outputs
  outputs = test->Evaluate(a, b);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 0);
}

TEST_F(WorkPieceTests, FixedInputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedInsMod>(3);

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate with a bool that will return both outputs
  obj->flag = true;
  std::vector<boost::any> outputs = test->Evaluate(s, a, obj);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate with a bool that will return one output
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), a);
}

TEST_F(WorkPieceTests, FixedOutputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedOutsMod>(2);

  // make sure the number of inputs and outputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // evaluate with no inputs
  std::vector<boost::any> outputs = test->Evaluate();

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare("")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);

  // evaluate with one inputs
  outputs = test->Evaluate(s);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, FixedInputOutputNum) {
  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(3, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);
}

TEST_F(WorkPieceTests, FixedInputTypes) {
  // the input types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInsMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), a);
}

TEST_F(WorkPieceTests, FixedOutputTypes) {
  // the output types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedOutsMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);
}

TEST_F(WorkPieceTests, SomeFixedInputTypes) {
  // the input types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[2] = typeid(std::shared_ptr<AnObject>).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInsMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj, 9.0, 'c');

  // make sure the outputs are what we expect
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), a);
}

TEST_F(WorkPieceTests, SomeFixedOutputTypes) {
  // the output types
  std::map<unsigned int, std::string> types;
  types[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedOutsMod>(types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  // evaluate with no inputs
  std::vector<boost::any> outputs = test->Evaluate(s);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  outputs = test->Evaluate(s, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), b);
}

TEST_F(WorkPieceTests, SomeFixedInputTypesFixedInputNum) {
    // the input types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[2] = typeid(std::shared_ptr<AnObject>).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInsMod>(types, 3);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a*b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[0]), a);
}

TEST_F(WorkPieceTests, SomeFixedInputTypesFixedOutputNum) {
  // the input types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[2] = typeid(std::shared_ptr<AnObject>).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(types, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj, 9.0, 'c');

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  outputs = test->Evaluate(s);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedOutputTypesFixedInputNumber) {
  // the output types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(3, types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedOutputTypesFixedOutputNumber) {
  // the output types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedOutsMod>(types, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, obj, a);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), (double)b);
}

TEST_F(WorkPieceTests, FixedInputTypesFixedOutputNum) {
  // the input types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(types, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 3 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, FixedOutputTypesFixedInputNum) {
  // the output types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(3, types);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 3 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedInputTypesFixedInputOutputNum) {
  // the input types
  std::map<unsigned int, std::string> types;
  types[0] = typeid(std::string).name();
  types[2] = typeid(std::shared_ptr<AnObject>).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(types, 3, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedOutputTypesFixedInputOutputNumber) {
  // the output types
  std::map<unsigned int, std::string> types;
  types[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(3, types, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, FixedTypes) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedInputTypesFixedOutputTypes) {
  // the input types
  std::map<unsigned int, std::string> inTypes;
  inTypes[0] = typeid(std::string).name();
  inTypes[1] = typeid(double).name();

  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj, 'v');

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj, 1, 9.0);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, SomeFixedInputTypesFixedInputNumFixedOutputTypes) {
  // the input types
  std::map<unsigned int, std::string> inTypes;
  inTypes[0] = typeid(std::string).name();
  inTypes[1] = typeid(double).name();

  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, 3, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);
}

TEST_F(WorkPieceTests, FixedInputTypesSomeFixedOutputTypes) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // the output types
  std::map<unsigned int, std::string> outTypes;
  outTypes[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(b);

  // evaluate the WorkPiece
  obj->flag = true;
  auto outputs = test->Evaluate(s, a, obj);

  // make sure we get 2 outputs
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), a);

  // evaluate the WorkPiece
  obj->flag = false;
  outputs = test->Evaluate(s, a, obj);

  // make sure we get 1 output
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs[0]).compare(s)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs[1]), 1.0);
}

TEST_F(WorkPieceTests, FixedInputTypesSomeFixedOutputTypesFixedOutputNum) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});

  // the output types
  std::map<unsigned int, std::string> outTypes;
  outTypes[1] = typeid(double).name();

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, outTypes, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);
}

TEST_F(WorkPieceTests, SomeFixedInOuts) {
  // the input types
  std::map<unsigned int, std::string> inTypes;
  inTypes[0] = typeid(std::string).name();
  inTypes[1] = typeid(double).name();

  // the output types
  std::map<unsigned int, std::string> outTypes;
  outTypes[1] = typeid(double).name();;

  // create the test WorkPiece
  auto test = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, -1);

  test = std::make_shared<FixedInOutMod>(inTypes, 3, outTypes);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, -1);

  test = std::make_shared<FixedInOutMod>(inTypes, outTypes, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, -1);
  EXPECT_EQ(test->numOutputs, 2);

  test = std::make_shared<FixedInOutMod>(inTypes, 3, outTypes, 2);

  // make sure the number of inputs matches what we expect
  EXPECT_EQ(test->numInputs, 3);
  EXPECT_EQ(test->numOutputs, 2);
}
