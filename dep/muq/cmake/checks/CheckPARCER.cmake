set(CMAKE_REQUIRED_LIBRARIES ${PARCER_LIBRARIES} pthread)
set(CMAKE_REQUIRED_INCLUDES ${PARCER_INCLUDE_DIRS})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <mpi.h>
  #include <parcer/Queue.h>

  using namespace parcer;

  class Work {
  public:
   inline Work() {};

   double Evaluate(double x) { return x*x; };
  };

  int main(int argc, char* argv[]){
    // Get the rank of the current process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::shared_ptr<Work> work = std::make_shared<Work>();
    auto comm = std::make_shared<Communicator>();
  }
  "
  PARCER_COMPILES)

if (NOT PARCER_COMPILES)
  set(PARCER_TEST_FAIL 1)
else ()
  set(PARCER_TEST_FAIL 0)
endif()
