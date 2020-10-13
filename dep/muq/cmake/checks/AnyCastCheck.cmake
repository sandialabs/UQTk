set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")

# try to compile a simple program to make sure the c++11 auto feature works
CHECK_CXX_SOURCE_COMPILES(
  "
class Tester {
public:
  Tester(double const& tempIn) : temp(tempIn){};

  template<typename T>
  operator T const&(){return temp;};

  double const& temp;
};

int main(){
  double a1 = 1.1;
  double a2 = Tester(a1);
  return 0;
}
  "
  ANYCAST_COMPILES)

if(ANYCAST_COMPILES)
  set(MUQ_ANYCAST_COMPILES 1  CACHE STRING "Result of ANYCAST_COMPILES test.")
else()
  set(MUQ_ANYCAST_COMPILES 0  CACHE STRING "Result of ANYCAST_COMPILES test.")  
endif()
