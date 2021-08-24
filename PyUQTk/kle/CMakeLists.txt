# Using Pybind to build the python module

find_package(pybind11)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
include_directories(../../cpp/lib/kle) # kle directory
pybind11_add_module(_kle PyKle.cpp)
TARGET_LINK_LIBRARIES(_kle PUBLIC uqtk)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/)

INSTALL(TARGETS _kle DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/kle/_kle.cpython-37m-darwin.so
	      DESTINATION PyUQTk/)