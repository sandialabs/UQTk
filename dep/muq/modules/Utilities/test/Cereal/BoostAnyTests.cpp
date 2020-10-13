
#include "MUQ/config.h"

#if MUQ_HAS_PARCER

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "MUQ/Utilities/Cereal/BoostAnySerializer.h"
#include <cereal/archives/binary.hpp>



TEST(Cerealize, BoostAnyDouble) {

  double temp = 10.0;
  boost::any startAny = temp;

  {
    std::ofstream out("eigen.cereal", std::ios::binary);
    cereal::BinaryOutputArchive archive_o(out);
    archive_o(startAny);
  }

  boost::any loadAny;
  {
    std::ifstream in("eigen.cereal", std::ios::binary);
    cereal::BinaryInputArchive archive_i(in);
    archive_i(loadAny);
  }

  EXPECT_DOUBLE_EQ(boost::any_cast<double>(startAny), boost::any_cast<double>(loadAny));
}

#endif
