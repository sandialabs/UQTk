#include <gtest/gtest.h>

#include "MUQ/Utilities/HDF5/PathTools.h"

TEST(Paths, GetParent)
{
    std::string base = "/my/base/path";
    std::string parent = muq::Utilities::GetParentPath(base);
    
    EXPECT_EQ("/my/base", parent);

    base += '/';
    parent = muq::Utilities::GetParentPath(base);
    EXPECT_EQ("/my/base", parent);
}
