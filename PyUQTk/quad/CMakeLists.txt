# Enable support for modern IN_LIST operator
if (POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

include_directories(../../cpp/lib/quad) #quad directory
include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/include) # include directory
pybind11_add_module(_quad pyquad.cpp)
TARGET_LINK_LIBRARIES(_quad PUBLIC uqtk blas lapack depuqtk)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/quad)

INSTALL(TARGETS _quad DESTINATION PyUQTk/)
