project (UQTk)

add_executable (GradHessianTest main.cpp)

target_link_libraries (GradHessianTest uqtk )

target_link_libraries (GradHessianTest depdsfmt )
target_link_libraries (GradHessianTest depcvode )
target_link_libraries (GradHessianTest depnvec  )
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
include_directories(../../../dep/cvode-2.7.0/include)

include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS GradHessianTest DESTINATION bin/tests/)

add_test(GradHessianTest GradHessianTest)
