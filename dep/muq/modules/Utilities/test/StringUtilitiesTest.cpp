
#include "gtest/gtest.h"

#include <iostream>
#include <fstream>

#include "MUQ/Utilities/StringUtilities.h"

using namespace std;
using namespace muq::Utilities;
using namespace muq::Utilities::StringUtilities;

TEST(Utilities, StringStrip)
{
    string truth = "String1";

    string base = "  String1";
    EXPECT_EQ(truth,Strip(base));

    base = "\t String1";
    EXPECT_EQ(truth,Strip(base));

    base = "\t\n\t   String1";
    EXPECT_EQ(truth,Strip(base));

    base = "String1  ";
    EXPECT_EQ(truth,Strip(base));

    base = "String1  \t";
    EXPECT_EQ(truth,Strip(base));

    base = "String1  \t\n";
    EXPECT_EQ(truth,Strip(base));

    base = "  \t\t String1  \t\n";
    EXPECT_EQ(truth,Strip(base));
}


TEST(Utilities, StringSplit)
{
    std::vector<std::string> truth = {"String1", "String2", "String3"};

    std::string base = "  String1, String2  ,  String3";

    auto parts = Split(base,',');
    EXPECT_EQ(3,parts.size());
    for(int i=0; i<parts.size(); ++i)
        EXPECT_EQ(truth.at(i), parts.at(i));

    base = "  String1 - String2  -  String3";

    parts = Split(base,'-');
    EXPECT_EQ(3,parts.size());
    for(int i=0; i<parts.size(); ++i)
        EXPECT_EQ(truth.at(i), parts.at(i));

}
