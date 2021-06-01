# Using Pybind to build the python module

find_package(pybind11)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/tools) #tools directory
include_directories(../../cpp/lib/include) #include directory
include_directories(../../dep/dsfmt) #dsfmt directory
include_directories(../../dep/figtree) #figtree directory


pybind11_add_module(_tools pyTools.cpp)
TARGET_LINK_LIBRARIES(_tools PUBLIC uqtk depuqtk)

SET(copy_FILES
  __init__.py
  )

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/tools)

INSTALL(TARGETS _tools DESTINATION PyUQTk/tools)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/tools/_tools.cpython-37m-darwin.so
	      DESTINATION PyUQTk/tools)
