set(CMAKE_REQUIRED_LIBRARIES ${GTEST_LIBRARIES} pthread)
set(CMAKE_REQUIRED_INCLUDES ${GTEST_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #define GTEST_USE_OWN_TR1_TUPLE 1
  #include <gtest/gtest.h>
  TEST(CMAKE_TEST, Dummy){
  EXPECT_EQ(1,1);
  EXPECT_NEAR(1.0,2.0,3.0);
  }
  int main(int argc, char* argv[]){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  }
  "
  GTEST_COMPILES)

if (NOT GTEST_COMPILES)
  set(GTEST_TEST_FAIL 1)
else ()
  set(GTEST_TEST_FAIL 0)
endif()
