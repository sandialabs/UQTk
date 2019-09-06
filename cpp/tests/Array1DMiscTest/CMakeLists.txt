project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
link_directories (${INSTALL_LIB_DIR})
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/cvode)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/nvector/serial)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/sunlinsol/dense)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/sunmatrix/dense)
add_executable (Array1DMiscTest main.cpp)

target_link_libraries (Array1DMiscTest uqtk )

target_link_libraries (Array1DMiscTest depdsfmt )
target_link_libraries (Array1DMiscTest sundials_cvode)
target_link_libraries (Array1DMiscTest sundials_nvecserial)
target_link_libraries (Array1DMiscTest sundials_sunlinsoldense)
target_link_libraries (Array1DMiscTest sundials_sunmatrixdense)
target_link_libraries (Array1DMiscTest depslatec)
target_link_libraries (Array1DMiscTest deplapack)
target_link_libraries (Array1DMiscTest depblas  )
target_link_libraries (Array1DMiscTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (Array1DMiscTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (Array1DMiscTest ifcore ifport)
  else()
    target_link_libraries (Array1DMiscTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (Array1DMiscTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (Array1DMiscTest gfortran stdc++)
  else()
    target_link_libraries (Array1DMiscTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/mcmc  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lapack)
include_directories(../../../dep/blas)
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
#include_directories(../../../dep/cvode-2.7.0/include)

#include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS Array1DMiscTest DESTINATION bin/tests/)

add_test(Array1DMiscTest Array1DMiscTest)
