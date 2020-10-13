
#include "MUQ/config.h"


#if MUQ_HAS_PARCER

#include <gtest/gtest.h>
#include <parcer/Eigen.h>
#include <cereal/archives/binary.hpp>

#include <iostream>
#include <fstream>

// adapted from: https://stackoverflow.com/questions/22884216/serializing-eigenmatrix-using-cereal-library

TEST(CerealizeEigen, RandomMatrix) {
  const Eigen::MatrixXd test = Eigen::MatrixXd::Random(10, 3);

  const std::string tempName = "eigen.cereal";

  {
    std::ofstream out(tempName, std::ios::binary);
    cereal::BinaryOutputArchive archive_o(out);
    archive_o(test);
  }

  Eigen::MatrixXd test_loaded;

  {
    std::ifstream in(tempName, std::ios::binary);
    cereal::BinaryInputArchive archive_i(in);
    archive_i(test_loaded);
  }

  EXPECT_DOUBLE_EQ((test-test_loaded).norm(), 0.0);

  std::remove(tempName.c_str());
}

#endif
