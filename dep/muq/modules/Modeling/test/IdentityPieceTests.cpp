#include "gtest/gtest.h"

#include "MUQ/Modeling/IdentityPiece.h"

using namespace muq::Modeling;

TEST(IdentityPiece, Unfixed) {
  // unknown inputs/outputs
  auto identity = std::make_shared<IdentityPiece>();

  // evaulate with one input
  auto result = identity->Evaluate(1.0);

  EXPECT_EQ(result.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(0)), 1.0);

  // evaluate with two inputs
  const std::string str = "string";
  result = identity->Evaluate(2.0, str);

  EXPECT_EQ(result.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(0)), 2.0);
  EXPECT_TRUE(boost::any_cast<std::string>(result.at(1)).compare(str)==0);
}

TEST(IdentityPiece, FixedNum) {
  // known input/output number
  auto identity = std::make_shared<IdentityPiece>(1);

  // evaulate with one input
  auto result = identity->Evaluate(1.0);

  EXPECT_EQ(result.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(0)), 1.0);

  // evaluate with two inputs
  const std::string str = "string";
  result = identity->Evaluate(str);

  EXPECT_EQ(result.size(), 1);
  EXPECT_TRUE(boost::any_cast<std::string>(result.at(0)).compare(str)==0);
}

TEST(IdentityPiece, FixedType) {
  // the types
  std::vector<std::string> types({typeid(std::string).name(), typeid(double).name()});

  // known input/output number
  auto identity = std::make_shared<IdentityPiece>(types);

  // evaulate 
  const std::string str = "string";
  auto result = identity->Evaluate(str, 1.0);

  EXPECT_EQ(result.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(result.at(0)).compare(str)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(1)), 1.0);
}

TEST(IdentityPiece, SomeFixedTypes) {
  // the types
  std::map<unsigned int, std::string> types;
  types[1] = typeid(double).name();

  // known input/output number
  auto identity = std::make_shared<IdentityPiece>(types);

  // evaulate 
  auto result = identity->Evaluate(1.0);

  EXPECT_EQ(result.size(), 1);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(0)), 1.0);

  // evaulate 
  const std::string str = "string";
  result = identity->Evaluate(str, 1.0);

  EXPECT_EQ(result.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(result.at(0)).compare(str)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(1)), 1.0);
}

TEST(IdentityPiece, SomeFixedTypesFixedNum) {
  // the types
  std::map<unsigned int, std::string> types;
  types[1] = typeid(double).name();

  // known input/output number
  auto identity = std::make_shared<IdentityPiece>(types, 2);

  // evaulate 
  auto result = identity->Evaluate(1.0, 2.0);

  EXPECT_EQ(result.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(0)), 1.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(1)), 2.0);

  // evaulate 
  const std::string str = "string";
  result = identity->Evaluate(str, 1.0);

  EXPECT_EQ(result.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(result.at(0)).compare(str)==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result.at(1)), 1.0);
}
