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
  UTILITIES_CORE_PYTHON
    "Python wrappers for the basic building blocks of the Utilities library."
    "pymuqUtilities"
    "UTILITIES_CORE"
    "PYTHON"
    ""
    Module.cpp
    #LinearAlgebraWrapper.cpp
    MultiIndicesWrapper.cpp
    GeneralUtilitiesWrapper.cpp
)
