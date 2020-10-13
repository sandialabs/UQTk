#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"

using namespace std;

using namespace muq::Utilities;

/*
  MultiIndexLimiter.MaxLimiter
  ---------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  limiter and false value for index outside the limiter.

  Test:
  Create MultiIndex at (0,0) and max order limiter of 2. Check to make sure the
  MultiIndex is feasible within the limiter.

    -----------                -----------
  2 |         |              2 |         |
    |         |                |         |
  1 |         |     --->     1 |         |
    |         |                |         |
  0 | x       |              0 | x       |
    -----------                -----------
      0   1   2                  0   1   2
*/
TEST(Utilities_MultiIndicies, MaxLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = make_shared<MaxOrderLimiter>(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,3);
  multi->SetValue(1,3);

  // Check the result of IsFeasible().
  feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}

/*
  MultiIndexLimiter.ValidTotalLimiter
  -----------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  total order limiter and false outside.

  Test:
  Create MultiIndex at (0,0) and total order limiter of 2. Check to make sure the
  MultiIndex is feasible within the limiter.


  2 |\                       2 |\
    |  \                       |  \
  1 |    \          --->     1 |    \
    |      \                   |      \
  0 | x      \               0 | x      \
    -----------                -----------
      0   1   2                  0   1   2
*/
TEST(Utilities_MultiIndicies, TotalLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = make_shared<TotalOrderLimiter>(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,2);

  // Check the result of IsFeasible().
  feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}


/*
  MultiIndexLimiter.ValidAndLimiter
  ---------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  and limiter.

  Test:
  Create MultiIndex at (0,0) with total order limiter of 2 and max order limiter
  of 2. Check to make sure the MultiIndex is feasible within both limiters.

    ------------                -----------
  2 |\         |              2 |\         |
    |  \       |                |  \       |
  1 |    \     |     --->     1 |    \     |
    |      \   |                |      \   |
  0 | x      \ |              0 | x      \ |
    -----------                 ------------
      0   1   2                  0   1   2
*/
TEST(Utilities_MultiIndicies, AndLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> totalLimiter = make_shared<TotalOrderLimiter>(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = make_shared<MaxOrderLimiter>(2);

  // Combine into AndLimiter.
  shared_ptr<MultiIndexLimiter> andLimiter = make_shared<AndLimiter>(totalLimiter, maxLimiter);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = andLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,2);

  // Check the result of IsFeasible().
  feasible = andLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,3);
  multi->SetValue(1,2);

  // Check the result of IsFeasible().
  feasible = andLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}


/*
  MultiIndexLimiter.ValidOrLimiter1
  ---------------------------------

  Purpose:
  Make sure the OrLimiter class returns True for case where both constituent
  limiters are True.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  OrLimiter with both. Create a MultiIndex at (1,1) and make sure it is
  feasible - should be True for both constituent limiters.
*/
TEST(Utilities_MultiIndicies, OrLimiter)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = make_shared<TotalOrderLimiter>(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = make_shared<MaxOrderLimiter>(2);

  // Combine into OrLimiter.
  shared_ptr<MultiIndexLimiter> orLimiter = make_shared<OrLimiter>(totalLimiter, maxLimiter);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,1);
  multi->SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = orLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,3);

  // Check the result of IsFeasible().
  feasible = orLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);


  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,1);

  // Check the result of IsFeasible().
  feasible = orLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}


/*
  MultiIndexLimiter.ValidXorLimiter
  ---------------------------------

  Purpose:
  Make sure the XorLimiter class returns False for case where one constituent
  limiter is True, and the other is False.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  XorLimiter with both. Create a MultiIndex at (2,1) and make sure it is
  feasible - should be True for MaxOrderLimiter, False for TotalOrderLimiter ->
  XOR = True.
*/
TEST(Utilities_MultiIndicies, XorLimiter)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = make_shared<TotalOrderLimiter>(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = make_shared<MaxOrderLimiter>(2);

  // Combine into XorLimiter.
  shared_ptr<MultiIndexLimiter> xorLimiter = make_shared<XorLimiter>(totalLimiter, maxLimiter);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = xorLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);


  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,1);
  multi->SetValue(1,1);

  // Check the result of IsFeasible().
  feasible = xorLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);

  multi = make_shared<MultiIndex>(2);
  multi->SetValue(0,2);
  multi->SetValue(1,3);

  // Check the result of IsFeasible().
  feasible = xorLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
// /*
//   MultiIndexLimiter.ValidGeneralLimiter
//   -------------------------------------

//   Purpose:
//   Make sure the GeneralLimiter class returns True value for index within the
//   GeneralLimiter.

//   Test:
//   Create MultiIndexSet with max order 1 and dimension 2, and input it to a
//   GeneralLimiter. Create MultiIndex at (0,0) and max order of 2. Check to make
//   sure the MultiIndex is feasible within the GeneralLimiter.

//   2 |                       2 |
//     |                         |
//   1 | x   x         --->    1 |
//     |                         |
//   0 | x   x                 0 | x
//     -----------               -----------
//       0   1   2                 0   1   2
// */

// TEST(Utilities_MultiIndicies, ValidGeneralLimiter)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // GeneralLimiter.
//   shared_ptr<MultiIndexLimiter> limiter = GeneralLimiter(indexFamily);

//   // MultiIndex for testing.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

//   // Check the result of IsFeasible().
//   bool feasible = limiter->IsFeasible(multi);
//   EXPECT_TRUE(feasible);
// }

// /*
//   MultiIndexLimiter.InvalidGeneralLimiter
//   ---------------------------------------

//   Purpose:
//   Make sure the GeneralLimiter class returns False value for index outside the
//   GeneralLimiter.

//   Test:
//   Create MultiIndexSet with max order 1 and dimension 2, and input it to a
//   GeneralLimiter. Create MultiIndex at (2,1) and max order of 2. Check to make
//   sure the MultiIndex is NOT feasible within the GeneralLimiter.

//   2 |                       2 |
//     |                         |
//   1 | x   x         --->    1 |         o
//     |                         |
//   0 | x   x                 0 |
//     -----------               -----------
//       0   1   2                 0   1   2
// */
// TEST(Utilities_MultiIndicies, InvalidGeneralLimiter)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // GeneralLimiter.
//   shared_ptr<MultiIndexLimiter> limiter = GeneralLimiter(indexFamily);

//   // MultiIndex for testing.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,2);
//   multi.SetValue(1,1);

//   // Check the result of IsFeasible().
//   bool feasible = limiter->IsFeasible(multi);
//   EXPECT_FALSE(feasible);
// }

/*
  MultiIndexLimiter.ValidDimensionLimiter
  ---------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns True value for index within the
  DimensionLimiter.

  Test:
  Create DimensionLimiter with lower dimension of 1 and length of 1, and create
  a MultiIndex of [1 0 0 1]. Check to make sure the MultiIndex is feasible
  within the DimensionLimiter.

  DimensionLimiter(1,1) of [a0  a1  a2  a3] =

*/
TEST(Utilities_MultiIndicies, ValidDimensionLimiter)
{
  // DimensionLimiter.
  shared_ptr<MultiIndexLimiter> limiter = make_shared<DimensionLimiter>(1,1);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(4);
  multi->SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);

  multi = make_shared<MultiIndex>(4);
  multi->SetValue(1,1);
  multi->SetValue(3,1);

  // Check the result of IsFeasible().
  feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
