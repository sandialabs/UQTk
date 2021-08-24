# Using Pybind to build the python module for PCE

find_package(pybind11)

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


if(BUILD_SUNDIALS)
	TARGET_LINK_LIBRARIES(_pce PUBLIC uqtk depuqtk blas lapack gfortran ${PROJECT_BINARY_DIR}/../dep/sundials/src/cvode/libsundials_cvode.a ${PROJECT_BINARY_DIR}/../dep/sundials/src/nvector/serial/libsundials_nvecserial.a ${PROJECT_BINARY_DIR}/../dep/sundials/src/sunlinsol/dense/libsundials_sunlinsoldense.a ${PROJECT_BINARY_DIR}/../dep/sundials/src/sunmatrix/dense/libsundials_sunmatrixdense.a)
else()
	TARGET_LINK_LIBRARIES(_pce PUBLIC uqtk depuqtk blas lapack gfortran ${CMAKE_SUNDIALS_DIR}/lib/libsundials_cvode.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_nvecserial.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_sunlinsoldense.a ${CMAKE_SUNDIALS_DIR}/lib/libsundials_sunmatrixdense.a)
endif()

#SET(copy_FILES
#  __init__.py
#  )

#INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/pce)

INSTALL(TARGETS _pce DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pce/_pce.cpython-37m-darwin.so DESTINATION PyUQTk/)
