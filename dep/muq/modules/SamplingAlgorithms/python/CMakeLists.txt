# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    SAMPLINGALGORITHMS_CORE_PYTHON
    "Python wrappers for the basic building blocks of the SamplingAlgorithms library."
    pymuqSamplingAlgorithms
    "SAMPLING_ALGORITHM"
    "PYTHON"
    ""
    Module.cpp
    KernelWrapper.cpp
    ProposalWrapper.cpp
    SampleWrapper.cpp
    MCMCWrapper.cpp
    ProblemWrapper.cpp
)
