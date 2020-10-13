#include <gtest/gtest.h>

#include <memory>

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Modeling;

TEST(SundialslgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // if muq was compiled with Sundials, check for Sundials types
#if MUQ_HAS_SUNDIALS==1
  N_Vector vec = N_VNew_Serial(8);
  for( unsigned int i=0; i<8; ++i ) { NV_Ith_S(vec, i) = 3.0*i; }

  EXPECT_EQ(alg->Size(vec), 8);

  N_VDestroy(vec);
#endif
}

TEST(SundialsAlgebraTests, AccessElement) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // if muq was compiled with Sundials, check for Sundials types
#if MUQ_HAS_SUNDIALS==1
  N_Vector vec = N_VNew_Serial(8);
  for( unsigned int i=0; i<8; ++i ) { NV_Ith_S(vec, i) = 3.0*i; }

  for( unsigned int i=0; i<8; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(vec, i)), NV_Ith_S(vec, i));
  }

  N_VDestroy(vec);
#endif
}
