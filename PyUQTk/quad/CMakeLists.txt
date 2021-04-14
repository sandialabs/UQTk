# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

IF(NOT pybind11_FOUND)
  add_subdirectory (../../dep/pybind11)
endif()

include_directories(../../cpp/lib/quad) #quad directory
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_quad pyquad.cpp)
TARGET_LINK_LIBRARIES(_quad PUBLIC uqtk blas lapack depuqtk)

INSTALL(TARGETS _quad DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/quad/_quad.cpython-37m-darwin.so DESTINATION PyUQTk)
