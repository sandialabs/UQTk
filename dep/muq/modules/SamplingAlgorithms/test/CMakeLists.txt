CreateTestGroup(
    SAMPLING_ALGORITHM

    SampleCollectionTests.cpp
    MarkovChainTests.cpp
    MIMCMCTests.cpp
    MLMCMCTests.cpp

    GMHKernelTests.cpp

    #MonteCarloTests.cpp
    ImportanceSamplingTests.cpp
    MCMCTests.cpp
    CrankNicolsonTests.cpp
)

CreateTestGroup(
  SAMPLING_PROBLEMS

  SamplingProblemTests.cpp
)

CreateTestGroup(
  FANCY_SAMPLING_PROBLEMS

  ExpensiveSamplingProblemTests.cpp
)

CreateParallelTestGroup(
    SAMPLING_ALGORITHM

    SamplingStateParallelTests.cpp

    DistributedCollectionTests.cpp

    GMHKernelParallelTests.cpp

    SingleChainMCMCParallelTests.cpp

    ExpensiveSamplingProblemParallelTests.cpp

    ImportanceSamplingParallelTests.cpp
)
