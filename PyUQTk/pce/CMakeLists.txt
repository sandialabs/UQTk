# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
include_directories(../../cpp/lib/pce) #pce directory
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_pce pypce.cpp)

INSTALL(TARGETS _pce DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pce/pce.py DESTINATION PyUQTk)
