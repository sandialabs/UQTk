# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(pyuqtkarray pyuqtkarray.cpp)
TARGET_LINK_LIBRARIES(pyuqtkarray PUBLIC uqtk blas lapack)

INSTALL(TARGETS pyuqtkarray DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pyuqtkarray/pyuqtkarray.py
	      DESTINATION PyUQTk)
