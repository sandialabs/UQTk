# Using Pybind to build the python module

if(NOT pybind11_FOUND)
  include_directories(../../dep/pybind11)
else()
  find_package(pybind11)
endif()

include_directories(../../cpp/lib/quad) #quad directory
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_quad pyquad.cpp)
TARGET_LINK_LIBRARIES(_quad PUBLIC uqtk blas lapack depuqtk)

SET(copy_FILES
  __init__.py
  )

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/quad)

INSTALL(TARGETS _quad DESTINATION PyUQTk/quad)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/quad/_quad.cpython-37m-darwin.so DESTINATION PyUQTk/quad)
