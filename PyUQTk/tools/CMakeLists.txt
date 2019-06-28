enable_language(Fortran)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH}/../../Extras/lib/python/numpy/core/include)

#include source files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/) # array classes, array input output, and array tools
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/include/) # utilities like error handlers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/) # tools like multindex, etc. 

# include dependencies
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/blas/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/lapack/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/dsfmt/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/figtree/) # blas library headers
# INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../numpy/) # numpy headers

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(tools.i PROPERTIES CPLUSPLUS ON)

# compile swig with cpp extensions
SWIG_ADD_MODULE(
	tools python tools.i 
	# array tools needed to compile misc tools source files
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arrayio.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arraytools.cpp

	# source code for tools
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/combin.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/gq.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/minmax.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/multiindex.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/pcmaps.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/probability.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/rosenblatt.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/toolsf.f
)

# link python and 3rd party libraries, e.g., gfortran and blas
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  SWIG_LINK_LIBRARIES(tools uqtkarray deplapack depdsfmt depblas depfigtree depann gfortran ${PYTHON_LIBRARIES})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  SWIG_LINK_LIBRARIES(tools uqtkarray deplapack depdsfmt depblas depfigtree depann ifcore ${PYTHON_LIBRARIES})
endif()


INSTALL(TARGETS _tools DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/tools/tools.py DESTINATION PyUQTk)
