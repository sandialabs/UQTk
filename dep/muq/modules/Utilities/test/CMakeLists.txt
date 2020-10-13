CreateTestGroup(
    UTILITIES_HDF5

    HDF5/DatasetTests.cpp
    HDF5/PathTests.cpp
    HDF5/HDF5FileTest.cpp
)

CreateTestGroup(
    UTILITIES_MULTIINDEX

    MultiIndices/MultiIndexTests.cpp
    MultiIndices/MultiIndexSetTests.cpp
    MultiIndices/MultiIndexLimiterTests.cpp
)

CreateTestGroup(
    UTILITIES_CORE

    RandomGeneratorTest.cpp
    AnyCastTest.cpp
    StringUtilitiesTest.cpp

    Cereal/CerealizeEigenTests.cpp
    Cereal/BoostAnyTests.cpp
)
