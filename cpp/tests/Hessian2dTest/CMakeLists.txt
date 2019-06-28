project (UQTk)

add_executable (Hessian2dTest main.cpp)

target_link_libraries (Hessian2dTest uqtk )

target_link_libraries (Hessian2dTest depdsfmt )
target_link_libraries (Hessian2dTest depcvode )
target_link_libraries (Hessian2dTest depnvec  )
target_link_libraries (Hessian2dTest depslatec)
target_link_libraries (Hessian2dTest deplapack)
target_link_libraries (Hessian2dTest depblas  )
target_link_libraries (Hessian2dTest deplbfgs  )
target_link_libraries (Hessian2dTest depfigtree  )
target_link_libraries (Hessian2dTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (Hessian2dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (Hessian2dTest ifcore ifport)
  else()
    target_link_libraries (Hessian2dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (Hessian2dTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (Hessian2dTest gfortran stdc++)
  else()
    target_link_libraries (Hessian2dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS Hessian2dTest DESTINATION bin/tests/)

add_test(Hessian2dTest Hessian2dTest)
