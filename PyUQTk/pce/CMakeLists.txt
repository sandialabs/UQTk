# Using Pybind to build the python module for PCE

find_package(pybind11 REQUIRED)

#include source files
INCLUDE_DIRECTORIES(../../cpp/lib/array/) # array classes, array input output, and array tools
INCLUDE_DIRECTORIES(../../cpp/lib/include/) # utilities like error handlers
INCLUDE_DIRECTORIES(../../cpp/lib/tools/) # tools like multindex, etc.
INCLUDE_DIRECTORIES(../../cpp/lib/quad/) # quad class
INCLUDE_DIRECTORIES(../../cpp/lib/kle/) # kle class
INCLUDE_DIRECTORIES(../../cpp/lib/pce/) # PCSet and PCBasis classes

# include dependencies
INCLUDE_DIRECTORIES(../../dep/dsfmt/) # dsfmt
INCLUDE_DIRECTORIES(../../dep/figtree/) # figtree
INCLUDE_DIRECTORIES(../../dep/slatec/) # slatec headers
# cvode
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
	include_directories ("${PROJECT_BINARY_DIR}/../dep/sundials/include")
	include_directories (../../dep/sundials/include )
endif()


pybind11_add_module(_pce pypce.cpp)

INSTALL(TARGETS _pce DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pce/pce.py DESTINATION PyUQTk)
