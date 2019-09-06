FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH}/../../Extras/lib/python/numpy/core/include)

#include source files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/) # array classes, array input output, and array tools
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/include/) # utilities like error handlers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/) # tools like multindex, etc.
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/quad/) # quad class
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/kle/) # kle class
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/pce/) # PCSet and PCBasis classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/bcs/) # bcs

# include dependencies
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/blas/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/lapack/) # blas library headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/dsfmt/) # dsfmt
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/figtree/) # figtree
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/slatec/) # slatec headers
include_directories (../../dep/sundials/include )
if( NOT BUILD_SUNDIALS)
	include_directories ("${PROJECT_BINARY_DIR}/../dep/sundials/include")
	include_directories (${CMAKE_SUNDIALS_DIR}/include)
endif()
INCLUDE_DIRECTORIES(${CMAKE_SUNDIALS_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../numpy/) # numpy headers

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(bcs.i PROPERTIES CPLUSPLUS ON)

# compile swig with cpp extensions
SWIG_ADD_LIBRARY(
	bcs LANGUAGE python SOURCES bcs.i
	# array tools needed to compile misc tools source files
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arrayio.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/array/arraytools.cpp

	# source code for quad and kle class
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/quad/quad.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/kle/kle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/bcs/bcs.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/pce/PCSet.cpp

	# source code for tools
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/combin.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/gq.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/minmax.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/multiindex.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/pcmaps.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/probability.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tools/rosenblatt.cpp
)

# link python and 3rd party libraries, e.g., gfortran and blas
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  SWIG_LINK_LIBRARIES(bcs uqtkpce uqtktools uqtkquad uqtkarray  deplapack depblas depslatec depdsfmt depann depfigtree gfortran ${PYTHON_LIBRARIES})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  SWIG_LINK_LIBRARIES(bcs uqtkpce uqtktools uqtkquad uqtkarray  deplapack depblas depslatec depdsfmt depann depfigtree ifcore ${PYTHON_LIBRARIES})
endif()

if(BUILD_SUNDIALS)
	SWIG_LINK_LIBRARIES(bcs ${PROJECT_BINARY_DIR}/../dep/sundials/src/cvode/libsundials_cvode.a)
	SWIG_LINK_LIBRARIES(bcs ${PROJECT_BINARY_DIR}/../dep/sundials/src/nvector/serial/libsundials_nvecserial.a)
	SWIG_LINK_LIBRARIES(bcs ${PROJECT_BINARY_DIR}/../dep/sundials/src/sunlinsol/dense/libsundials_sunlinsoldense.a)
	SWIG_LINK_LIBRARIES(bcs ${PROJECT_BINARY_DIR}/../dep/sundials/src/sunmatrix/dense/libsundials_sunmatrixdense.a)
else()
	SWIG_LINK_LIBRARIES(bcs ${CMAKE_SUNDIALS_DIR}/lib/libsundials_cvode.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_nvecserial.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_sunlinsoldense.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_sunmatrixdense.a)
endif()

INSTALL(TARGETS _bcs DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/bcs/bcs.py DESTINATION PyUQTk)
