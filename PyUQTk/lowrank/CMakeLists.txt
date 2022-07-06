# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

include_directories(../../cpp/lib/lowrank) #quad directory
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_lowrank pylowrank.cpp)
TARGET_LINK_LIBRARIES(_lowrank PUBLIC uqtk blas lapack depuqtk)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/lowrank)

INSTALL(TARGETS _lowrank DESTINATION PyUQTk/)
