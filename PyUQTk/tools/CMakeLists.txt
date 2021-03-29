# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/tools) #tools directory
include_directories(../../cpp/lib/include) #include directory
include_directories(../../dep/dsfmt) #dsfmt directory
include_directories(../../dep/figtree) #figtree directory


pybind11_add_module(_tools pyTools.cpp)
TARGET_LINK_LIBRARIES(_tools PUBLIC uqtk depuqtk)

INSTALL(TARGETS _tools DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/tools/tools.py
	      DESTINATION PyUQTk)
