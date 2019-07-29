project (UQTk)

add_executable (BCS1dTest main.cpp)

target_link_libraries (BCS1dTest uqtk )

target_link_libraries (BCS1dTest depdsfmt )
target_link_libraries (BCS1dTest depcvode )
target_link_libraries (BCS1dTest depnvec  )
target_link_libraries (BCS1dTest depslatec)
target_link_libraries (BCS1dTest deplapack)
target_link_libraries (BCS1dTest depblas  )
target_link_libraries (BCS1dTest deplbfgs  )
target_link_libraries (BCS1dTest depfigtree  )
target_link_libraries (BCS1dTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (BCS1dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (BCS1dTest ifcore ifport)
  else()
    target_link_libraries (BCS1dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (BCS1dTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (BCS1dTest gfortran stdc++)
  else()
    target_link_libraries (BCS1dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce    )
include_directories(../../lib/bcs    )
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

INSTALL(TARGETS BCS1dTest DESTINATION bin/tests/)

add_test(BCS1dTest BCS1dTest)
