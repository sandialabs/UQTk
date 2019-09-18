FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH}/../../Extras/lib/python/numpy/core/include)

#include source files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tmcmc/) # tmcmc

# include dependencies
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../dep/dsfmt/) # dsfmt

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(tmcmc.i PROPERTIES CPLUSPLUS ON)

# compile swig with cpp extensions
SWIG_ADD_LIBRARY(
	tmcmc LANGUAGE python SOURCES tmcmc.i
	${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/lib/tmcmc/tmcmc.cpp
)

# link python and 3rd party libraries, e.g., gfortran and blas
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  SWIG_LINK_LIBRARIES(tmcmc depslatec depdsfmt gfortran ${PYTHON_LIBRARIES})
	SWIG_LINK_LIBRARIES(tmcmc m lapack blas)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  SWIG_LINK_LIBRARIES(tmcmc depslatec depdsfmt ifcore ifport ${PYTHON_LIBRARIES})
	SWIG_LINK_LIBRARIES(tmcmc m lapack blas)
endif()

INSTALL(TARGETS _tmcmc DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/tmcmc/tmcmc.py DESTINATION PyUQTk)
