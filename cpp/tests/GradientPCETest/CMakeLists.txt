project (UQTk)

add_executable (GradientPCETest main.cpp)

target_link_libraries (GradientPCETest uqtk )

target_link_libraries (GradientPCETest depdsfmt )
target_link_libraries (GradientPCETest depcvode )
target_link_libraries (GradientPCETest depnvec  )
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

INSTALL(TARGETS GradientPCETest DESTINATION bin/tests/)

add_test(GradientPCETest GradientPCETest)
