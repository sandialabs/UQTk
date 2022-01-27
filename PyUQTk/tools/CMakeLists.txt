# Using Pybind to build the python module

find_package(pybind11 REQUIRED)

include_directories(../../cpp/lib/array) #array directory
include_directories(../../cpp/lib/tools) #tools directory
include_directories(../../cpp/lib/include) #include directory
include_directories(../../dep/dsfmt) #dsfmt directory
include_directories(../../dep/figtree) #figtree directory

include_directories (${CMAKE_SUNDIALS_DIR}/include)
include_directories (${CMAKE_SUNDIALS_DIR}/include/cvodes)
include_directories (${CMAKE_SUNDIALS_DIR}/include/cvode)
if( BUILD_SUNDIALS)
  #include_directories (${PROJECT_BINARY_DIR}/../../dep/sundials/include)
  include_directories (${PROJECT_BINARY_DIR}/../dep/sundials/include)
  include_directories ($ENV{UQTK_SRC}/dep/sundials/include/cvodes)
  include_directories ($ENV{UQTK_SRC}/dep/sundials/include/cvode)
  include_directories (../../dep/sundials/include)
endif()


pybind11_add_module(_tools pyTools.cpp)
target_link_directories(_tools PUBLIC "${PROJECT_BINARY_DIR}/../dep/sundials/src/cvode/")
TARGET_LINK_LIBRARIES(_tools PUBLIC uqtk depuqtk lapack blas sundials_cvode)

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/tools)

INSTALL(TARGETS _tools DESTINATION PyUQTk/)
