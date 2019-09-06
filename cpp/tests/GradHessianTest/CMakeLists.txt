project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/cvode")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/nvector/serial")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunlinsol/dense")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunmatrix/dense")

add_executable (GradHessianTest main.cpp)

target_link_libraries (GradHessianTest uqtk )

target_link_libraries (GradHessianTest depdsfmt )
target_link_libraries (GradHessianTest sundials_cvode)
target_link_libraries (GradHessianTest sundials_nvecserial)
target_link_libraries (GradHessianTest sundials_sunlinsoldense)
target_link_libraries (GradHessianTest sundials_sunmatrixdense)
target_link_libraries (GradHessianTest depslatec)
target_link_libraries (GradHessianTest deplapack)
target_link_libraries (GradHessianTest depblas  )
target_link_libraries (GradHessianTest deplbfgs  )
target_link_libraries (GradHessianTest depfigtree  )
target_link_libraries (GradHessianTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (GradHessianTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (GradHessianTest ifcore ifport)
  else()
    target_link_libraries (GradHessianTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (GradHessianTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (GradHessianTest gfortran stdc++)
  else()
    target_link_libraries (GradHessianTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS GradHessianTest DESTINATION bin/tests/)

add_test(GradHessianTest GradHessianTest)
