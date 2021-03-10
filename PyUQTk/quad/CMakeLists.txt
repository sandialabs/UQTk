# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
include_directories(../../cpp/lib/quad) #quad directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_quad pyquad.cpp)

INSTALL(TARGETS _quad DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/quad/quad.py DESTINATION PyUQTk)
