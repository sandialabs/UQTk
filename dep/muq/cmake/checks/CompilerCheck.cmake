set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")

# try to compile a simple program to make sure the c++11 auto feature works
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <vector>
  int main(){
   std::vector<double> vec0(2);
   auto vec = vec0;
    return 0; 
   }
  "
  CPP11_AUTO_COMPILES)


  # try to compile a simple program to make sure the c++11 shared pointer feature works
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <memory>
    int main(){
     std::shared_ptr<double> vec0;
	 vec0 = std::make_shared<double>(2.0);
      return 0; 
     }
    "
    CPP11_SHAREDPTR_COMPILES)
	
	
	if(NOT CPP11_AUTO_COMPILES OR NOT CPP11_SHAREDPTR_COMPILES)
		message( FATAL_ERROR "Basic tests of c++11 cannot be compiled.  Please make sure your compiler supports all c++11 features.\n" )
	endif()
