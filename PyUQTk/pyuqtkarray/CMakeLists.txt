# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

IF(NOT pybind11_FOUND)
  add_subdirectory (../../dep/pybind11)
endif()

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(pyuqtkarray pyuqtkarray.cpp)
TARGET_LINK_LIBRARIES(pyuqtkarray PUBLIC uqtk blas lapack)

SET(copy_FILES
  __init__.py
  )

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/pyuqtkarray)

INSTALL(TARGETS pyuqtkarray DESTINATION PyUQTk/pyuqtkarray)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pyuqtkarray/pyuqtkarray.cpython-37m-darwin.so
	      DESTINATION PyUQTk/pyuqtkarray)
