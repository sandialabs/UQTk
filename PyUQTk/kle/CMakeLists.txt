# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
include_directories(../../cpp/lib/kle) # include directory
pybind11_add_module(_kle PyKle.cpp)
TARGET_LINK_LIBRARIES(_kle PUBLIC uqtk)

INSTALL(TARGETS _kle DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/kle/kle.py
	      DESTINATION PyUQTk)