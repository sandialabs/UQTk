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
    APPROXIMATION_CORE_PYTHON
    "Python wrappers for the basic building blocks of the Approximation library."
    "pymuqApproximationWrappers"
    "APPROXIMATION_GP;APPROXIMATION_GP_Kernels"
    "PYTHON"
    ""
    wrappers/Module.cpp
    wrappers/KernelWrapper.cpp
    wrappers/GaussianWrapper.cpp
    wrappers/PolynomialsWrapper.cpp
    wrappers/KLWrapper.cpp
    wrappers/QuadratureWrapper.cpp
    wrappers/PolynomialChaosWrapper.cpp
)
