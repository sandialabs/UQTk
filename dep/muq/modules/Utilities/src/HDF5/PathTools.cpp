#include "MUQ/Utilities/HDF5/PathTools.h"

#include <iostream>
#include <assert.h>

std::string muq::Utilities::GetParentPath(std::string const& base)
{
    int lastInd = base.length()-1;
    if(base[lastInd]=='/')
   	lastInd--;
       
    int pos = base.find_last_of('/',lastInd);

    return base.substr(0,pos);
}



std::pair<std::string, std::string> muq::Utilities::SplitString(std::string const& path)
{
  assert(path.length()>1);
  int i;
  for(i=1; i<path.length(); ++i)
  {
    if(path[i]=='/')
      break;
  }
  std::string nextNode      = path.substr(0,i);
  std::string remainingPath = path.substr(i,path.length()-i);
  
  return make_pair(nextNode,remainingPath);
}
