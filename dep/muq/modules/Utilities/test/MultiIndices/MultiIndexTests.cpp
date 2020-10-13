
#include "gtest/gtest.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

using namespace std;
using namespace muq::Utilities;


TEST(Utilities_MultiIndices, DefaultConstructor)
{
  MultiIndex multi(2);

  EXPECT_EQ(2,multi.GetLength());
  EXPECT_EQ(0,multi.GetValue(0));
  EXPECT_EQ(0,multi.GetValue(1));
  EXPECT_EQ(0,multi.Max());
  EXPECT_EQ(0,multi.Sum());
}

TEST(Utilities_MultiIndices, CreateVector)
{
  Eigen::RowVectorXi truth(6);
  truth << 1,0,3, 0, 0, 1;

  MultiIndex multi(truth);

  Eigen::RowVectorXi test = multi.GetVector();
  for(int i=0; i<truth.size(); ++i)
    EXPECT_EQ(truth(i), test(i));
}

TEST(Utilities_MultiIndices, RowConstructor)
{
  Eigen::RowVectorXi temp(3);
  temp(0) = 1;
  temp(1) = 5;
  temp(2) = 0;

  MultiIndex multi(temp);

  EXPECT_EQ(3,multi.GetLength());
  EXPECT_EQ(1,multi.GetValue(0));
  EXPECT_EQ(5,multi.GetValue(1));
  EXPECT_EQ(5,multi.Max());
  EXPECT_EQ(6,multi.Sum());
}

TEST(Utilities_MultiIndices, InitializerListConstructor)
{
  MultiIndex multi{1,5,3,0};

  EXPECT_EQ(4,multi.GetLength());
  EXPECT_EQ(1,multi.GetValue(0));
  EXPECT_EQ(5,multi.GetValue(1));
  EXPECT_EQ(3,multi.GetValue(2));
  EXPECT_EQ(0,multi.GetValue(3));

  EXPECT_EQ(5,multi.Max());
  EXPECT_EQ(9,multi.Sum());
}

TEST(Utilities_MultiIndices, SetValue)
{

  // Start off with [0, 7]
  MultiIndex multi(2);
  multi.SetValue(1,7);

  EXPECT_EQ(0, multi.GetValue(0));
  EXPECT_EQ(7, multi.GetValue(1));

  EXPECT_EQ(7, multi.Max());
  EXPECT_EQ(7, multi.Sum());
  EXPECT_EQ(2, multi.GetLength());

  // [2, 7]
  multi.SetValue(0,2);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(7, multi.GetValue(1));
  EXPECT_EQ(7, multi.Max());
  EXPECT_EQ(9, multi.Sum());

  // [2,8]
  multi.SetValue(1,8);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(8, multi.GetValue(1));
  EXPECT_EQ(8, multi.Max());
  EXPECT_EQ(10, multi.Sum());

  // [2,3]
  multi.SetValue(1,3);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(3, multi.GetValue(1));
  EXPECT_EQ(3, multi.Max());
  EXPECT_EQ(5, multi.Sum());

 // [2, 0]
  multi.SetValue(1,0);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(0, multi.GetValue(1));
  EXPECT_EQ(2, multi.Max());
  EXPECT_EQ(2, multi.Sum());
}

TEST(Utilities_MultiIndices, AdditionOperators)
{
  MultiIndex a {1, 2, 3, 4};
  MultiIndex b {5, 6, 7, 8};
  MultiIndex c = a+b;

  EXPECT_EQ(1,a.GetValue(0));
  EXPECT_EQ(2,a.GetValue(1));
  EXPECT_EQ(3,a.GetValue(2));
  EXPECT_EQ(4,a.GetValue(3));

  EXPECT_EQ(5,b.GetValue(0));
  EXPECT_EQ(6,b.GetValue(1));
  EXPECT_EQ(7,b.GetValue(2));
  EXPECT_EQ(8,b.GetValue(3));

  EXPECT_EQ(6,c.GetValue(0));
  EXPECT_EQ(8,c.GetValue(1));
  EXPECT_EQ(10,c.GetValue(2));
  EXPECT_EQ(12,c.GetValue(3));

  a += b;

  EXPECT_EQ(6,a.GetValue(0));
  EXPECT_EQ(8,a.GetValue(1));
  EXPECT_EQ(10,a.GetValue(2));
  EXPECT_EQ(12,a.GetValue(3));

  EXPECT_EQ(5,b.GetValue(0));
  EXPECT_EQ(6,b.GetValue(1));
  EXPECT_EQ(7,b.GetValue(2));
  EXPECT_EQ(8,b.GetValue(3));

  ++b;

  EXPECT_EQ(6,b.GetValue(0));
  EXPECT_EQ(7,b.GetValue(1));
  EXPECT_EQ(8,b.GetValue(2));
  EXPECT_EQ(9,b.GetValue(3));
}

TEST(Utilities_MultiIndices, SubtractionOperators)
{
  MultiIndex a {1, 2, 3, 4};
  MultiIndex b {5, 6, 7, 8};
  MultiIndex c = b-a;

  EXPECT_EQ(1,a.GetValue(0));
  EXPECT_EQ(2,a.GetValue(1));
  EXPECT_EQ(3,a.GetValue(2));
  EXPECT_EQ(4,a.GetValue(3));

  EXPECT_EQ(5,b.GetValue(0));
  EXPECT_EQ(6,b.GetValue(1));
  EXPECT_EQ(7,b.GetValue(2));
  EXPECT_EQ(8,b.GetValue(3));

  EXPECT_EQ(4,c.GetValue(0));
  EXPECT_EQ(4,c.GetValue(1));
  EXPECT_EQ(4,c.GetValue(2));
  EXPECT_EQ(4,c.GetValue(3));

  b -= a;

  EXPECT_EQ(1,a.GetValue(0));
  EXPECT_EQ(2,a.GetValue(1));
  EXPECT_EQ(3,a.GetValue(2));
  EXPECT_EQ(4,a.GetValue(3));

  EXPECT_EQ(4,b.GetValue(0));
  EXPECT_EQ(4,b.GetValue(1));
  EXPECT_EQ(4,b.GetValue(2));
  EXPECT_EQ(4,b.GetValue(3));

  --a;
  --a;

  EXPECT_EQ(0,a.GetValue(0));
  EXPECT_EQ(0,a.GetValue(1));
  EXPECT_EQ(1,a.GetValue(2));
  EXPECT_EQ(2,a.GetValue(3));
}
