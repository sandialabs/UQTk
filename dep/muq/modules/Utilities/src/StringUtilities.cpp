#include "MUQ/Utilities/StringUtilities.h"

using namespace std;

std::vector<std::string> muq::Utilities::StringUtilities::Split(std::string str, char delim)
{
    std::vector<std::string> output;

    int pos = 0;
    std::string strPart;
    // while we keep finding commas...
    while ((pos = str.find(delim)) != std::string::npos) {

        // extract the part before the comma
        strPart = str.substr(0, pos);

        // erase the part before the comma
        str.erase(0, pos + 1);

        // Strip whitespace and store the section name
        output.push_back(muq::Utilities::StringUtilities::Strip(strPart));
    }

    // make sure to store whatevers left
    output.push_back(muq::Utilities::StringUtilities::Strip(str));

    return output;
}


std::string muq::Utilities::StringUtilities::Strip(std::string str)
{
    // remove whitespace from the beginning
    while((str.front()==' ')||(str.front()=='\t')||(str.front()=='\n'))
        str.erase(0,1);

    // remove whitespace from the end
    while((str.back()==' ')||(str.back()=='\t')||(str.back()=='\n'))
        str.erase(str.size()-1,1);

    return str;
}
