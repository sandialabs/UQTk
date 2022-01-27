# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/bcs) #bcs directory
include_directories(../../cpp/lib/include) #bcs directory
pybind11_add_module(_bcs PyBCS.cpp)
TARGET_LINK_LIBRARIES(_bcs PUBLIC uqtk depuqtk lapack blas)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk)

SET(copy_FILES
  bcs_ext.py
  )

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bcs_ext.py
               ${CMAKE_CURRENT_BINARY_DIR}/bcs_ext.py COPYONLY)




INSTALL(TARGETS _bcs DESTINATION PyUQTk/)
#INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/bcs/_bcs.cpython-37m-darwin.so
#	      DESTINATION PyUQTk/bcs)
