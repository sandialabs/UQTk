project (UQTk)

add_executable (GradHessianTest main.cpp)

target_link_libraries (GradHessianTest uqtk )

target_link_libraries (GradHessianTest depdsfmt )
target_link_libraries (GradHessianTest sundials_cvode)
target_link_libraries (GradHessianTest sundials_nvecserial)
target_link_libraries (GradHessianTest sundials_sunlinsoldense)
target_link_libraries (GradHessianTest sundials_sunmatrixdense)
target_link_libraries (GradHessianTest depslatec)
target_link_libraries (GradHessianTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (GradHessianTest m blas ${BLAS_LIBRARIES})
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
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
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
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()
INSTALL(TARGETS GradHessianTest DESTINATION bin/tests/)

add_test(GradHessianTest GradHessianTest)
