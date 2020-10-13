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
    UTILITIES_HDF5
    "MUQ's wrapper for interacting with HDF5 binary files."
    muqUtilities
    ""
    "EIGEN3;HDF5"
    ""
    HDF5/Attributes.cpp
    HDF5/H5Object.cpp
    HDF5/HDF5File.cpp
    HDF5/HDF5Types.cpp
    HDF5/PathTools.cpp
    HDF5/BlockDataset.cpp
)

CreateCompileGroup(
    UTILITIES_CORE
    "Basic utilties used throughout MUQ."
    muqUtilities
    ""
    "EIGEN3"
    ""
    RandomGenerator.cpp
    StringUtilities.cpp
)

CreateCompileGroup(
    UTILITIES_MULTIINDEX
    "Classes for defining and generating multiindices and multiindex sets."
    muqUtilities
    ""
    "EIGEN3"
    ""
    MultiIndices/MultiIndex.cpp
    MultiIndices/MultiIndexLimiter.cpp
    MultiIndices/MultiIndexSet.cpp
    MultiIndices/MultiIndexFactory.cpp
)

#CreateCompileGroup(
#    UTILITIES_PYTHONHELPERS_PYTHON
#    "Classes and functions for python helper functions (not bindings)."
#    pymuqUtilities
#    ""
#    "PYTHON"
#    ""
#    PyDictConversion.cpp
#)
