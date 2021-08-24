# Using Pybind to build the python module

find_package(pybind11)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/bcs) #bcs directory
include_directories(../../cpp/lib/include) #bcs directory
pybind11_add_module(_bcs PyBCS.cpp)
TARGET_LINK_LIBRARIES(_bcs PUBLIC uqtk)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk)


INSTALL(TARGETS _bcs DESTINATION PyUQTk/)
#INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/bcs/_bcs.cpython-37m-darwin.so
#	      DESTINATION PyUQTk/bcs)