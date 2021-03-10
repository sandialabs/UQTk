FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR}/../../Extras/lib/python/numpy/core/include)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/) # array classes, array input output, and array tools
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/include/) # utilities like error handlers

# INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../numpy/) # numpy headers

# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
pybind11_add_module(pyuqtkarray pyuqtkarray/pyuqtkarray.cpp)
TARGET_LINK_LIBRARIES(pyuqtkarray PUBLIC uqtkarray)

INSTALL(TARGETS pyuqtkarray DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pyuqtkarray/pyuqtkarray.py
	      DESTINATION PyUQTk)
