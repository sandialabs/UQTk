#include <iostream>

// include the google testing header
#include <gtest/gtest.h>

#include "MUQ/config.h"

#if MUQ_HAS_MPI==0
#error
#endif

#include <mpi.h>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  const int ierr = MPI_Init(nullptr, nullptr);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // set the random number generator (important to make sure Eigen::****::Random returns different values on each processor)
  srand((unsigned int) time(0)+rank);
  for( unsigned int i=0; i<100; ++i ) { rand(); }

  const int res = RUN_ALL_TESTS();   
  
  MPI_Finalize();
  
  return res;
}
