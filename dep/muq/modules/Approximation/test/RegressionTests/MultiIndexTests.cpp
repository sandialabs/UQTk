#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/MultiIndex.h"

using namespace muq::Approximation;

TEST(MultiIndex, Basic) {
  const unsigned int dim = 3;
  const unsigned int maxOrder = 2;
  const unsigned int minOrder = 0;
  
  // create a multi index
  auto multi = std::make_shared<MultiIndex>(dim, maxOrder, minOrder);

  unsigned int cnt = 1+dim+dim*(dim+1)/2;
  EXPECT_EQ(multi->Size(), 1+dim+dim*(dim+1)/2);

  // evaluate the multi-index
  for( unsigned int i=0; i<multi->Size(); ++i ) {
    // get the ith multi-index
    const std::vector<boost::any>& result = multi->Evaluate(i);
    const Eigen::VectorXi& alpha = boost::any_cast<Eigen::VectorXi const&>(result[0]);

    EXPECT_EQ(alpha.size(), dim);
    EXPECT_TRUE(alpha.sum()<=maxOrder);
    EXPECT_TRUE(alpha.sum()>=minOrder);
  }
}
