
# we want to decide if adding the -stdlib=libc++ compiler flag is a good thing
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <cstddef>
	#include <locale>
	#include <exception>
	#include <ostream>
	#include <istream>
	#include <string>
	#include <cstring>
	#include <algorithm>
	#include <set>
	#include <memory>
    int main(){
     std::shared_ptr<double> vec0;
	 vec0 = std::make_shared<double>(2.0);
      return 0; 
     }
    "
    CPP11_STDLIB_COMPILES)
	
	if(CPP11_STDLIB_COMPILES)
		set(MUQ_USE_LIBC11 ON)
		set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -stdlib=libc++")
		message(STATUS "The compiler was found to work with libc++.  Adding -stdlib=libc++ to CMAKE_CXX_FLAGS.")
	else()
		set(MUQ_USE_LIBC11 OFF)
		message(STATUS "The compiler was NOT found to work with libc++.")
		set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
	endif()