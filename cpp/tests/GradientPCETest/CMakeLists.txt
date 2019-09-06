project (UQTk)

add_executable (GradientPCETest main.cpp)

target_link_libraries (GradientPCETest uqtk )

target_link_libraries (GradientPCETest depdsfmt )
target_link_libraries (GradientPCETest sundials_cvode)
target_link_libraries (GradientPCETest sundials_nvecserial)
target_link_libraries (GradientPCETest sundials_sunlinsoldense)
target_link_libraries (GradientPCETest sundials_sunmatrixdense)
target_link_libraries (GradientPCETest depslatec)
target_link_libraries (GradientPCETest deplapack)
target_link_libraries (GradientPCETest depblas  )
target_link_libraries (GradientPCETest deplbfgs  )
target_link_libraries (GradientPCETest depfigtree  )
target_link_libraries (GradientPCETest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (GradientPCETest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (GradientPCETest ifcore ifport)
  else()
    target_link_libraries (GradientPCETest ${IntelLibPath}/libifcore.a)
    target_link_libraries (GradientPCETest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (GradientPCETest gfortran stdc++)
  else()
    target_link_libraries (GradientPCETest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

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
include_directories (${CMAKE_SUNDIALS_DIR}/include)

include_directories (${INSTALL_LIB_DIR})
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")
include_directories(/Users/kjohnst/UQTk-build/dep/sundials/include)

INSTALL(TARGETS GradientPCETest DESTINATION bin/tests/)

add_test(GradientPCETest GradientPCETest)
