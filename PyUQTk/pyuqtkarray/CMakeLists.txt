# Using Pybind to build the python module

find_package(pybind11)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_uqtkarray pyuqtkarray.cpp)
TARGET_LINK_LIBRARIES(_uqtkarray PUBLIC uqtk blas lapack)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/pyuqtkarray)

INSTALL(TARGETS _uqtkarray DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pyuqtkarray/_uqtkarray.cpython-37m-darwin.so
	      DESTINATION PyUQTk/)
