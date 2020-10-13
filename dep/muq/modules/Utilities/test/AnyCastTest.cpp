#include "gtest/gtest.h"

#include "MUQ/Utilities/AnyHelpers.h"

using namespace muq::Utilities;

TEST(UtilitiesAnyHelpers, AnyCast)
{
  double a = 1.123;
  boost::any aAny = a;

  double b = AnyCast(aAny);
  EXPECT_EQ(a,b);

  double& bRef = AnyCast(aAny);
  double& bRefBoost = boost::any_cast<double&>(aAny);
  EXPECT_EQ(a,bRef);
  EXPECT_EQ(&bRefBoost, &bRef);

  double const& bConstRef = AnyCast(aAny);
  double const& bConstRefBoost = boost::any_cast<double const&>(aAny);
  EXPECT_EQ(a,bConstRef);
  EXPECT_EQ(&bConstRefBoost, &bConstRef);
}

TEST(UtilitiesAnyHelpers, AnyConstCast)
{
  const double a = 1.123;
  const boost::any aAny = a;

  double b = AnyConstCast(aAny);
  EXPECT_EQ(a,b);

  double const& bRef = AnyConstCast(aAny);
  double const& bRefBoost = boost::any_cast<double const&>(aAny);
  EXPECT_EQ(a,b);
  EXPECT_EQ(&bRefBoost, &bRef);
}
