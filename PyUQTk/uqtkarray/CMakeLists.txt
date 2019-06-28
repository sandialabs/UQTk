FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR}/../../Extras/lib/python/numpy/core/include)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/) # array classes, array input output, and array tools
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/include/) # utilities like error handlers

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/blas/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/lapack/) # blas library headers
# INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../numpy/) # numpy headers

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(uqtkarray.i PROPERTIES CPLUSPLUS ON)

# compile swig with cpp extensions
SWIG_ADD_MODULE(uqtkarray python uqtkarray.i 
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arrayio.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arraytools.cpp
)

# link python and 3rd party libraries, e.g., gfortran and blas
# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  SWIG_LINK_LIBRARIES(uqtkarray deplapack depblas gfortran ${PYTHON_LIBRARIES})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  SWIG_LINK_LIBRARIES(uqtkarray deplapack depblas ifcore ${PYTHON_LIBRARIES})
endif()


INSTALL(TARGETS _uqtkarray DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/uqtkarray/uqtkarray.py
	      DESTINATION PyUQTk)
