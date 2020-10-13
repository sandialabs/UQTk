# make sure that the boost graph library is available
set(CMAKE_REQUIRED_LIBRARIES ${BOOST_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
"
#include <boost/graph/adjacency_list.hpp>
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS, int, int> Graph;
int main(){
Graph temp;
return 0;
}
"
BOOST_GRAPH_COMPILES)


CHECK_CXX_SOURCE_COMPILES(
"
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/filesystem.hpp>
using namespace std;
using namespace boost::filesystem;
int main(int argc, char* argv[])
{
path p (argv[1]);
directory_iterator temp(p);
directory_iterator temp2(temp);
return 0;
}
"
BOOST_DIRECTORY_ITERATOR_COMPILES)

if(NOT BOOST_GRAPH_COMPILES OR NOT BOOST_DIRECTORY_ITERATOR_COMPILES)
	set(BOOST_TEST_FAIL 1)
else()
	set(BOOST_TEST_FAIL 0)
endif()
