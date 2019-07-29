project (UQTk)

add_executable (ArraySortTest main.cpp)

target_link_libraries (ArraySortTest uqtk  )

target_link_libraries (ArraySortTest depdsfmt )
target_link_libraries (ArraySortTest depcvode )
target_link_libraries (ArraySortTest depnvec  )
target_link_libraries (ArraySortTest depslatec)
target_link_libraries (ArraySortTest deplapack)
target_link_libraries (ArraySortTest depblas  )
target_link_libraries (ArraySortTest deplbfgs  )
target_link_libraries (ArraySortTest depfigtree )
target_link_libraries (ArraySortTest  depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (ArraySortTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (ArraySortTest ifcore)
  else()
    target_link_libraries (ArraySortTest ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (ArraySortTest gfortran stdc++)
  else()
    target_link_libraries (ArraySortTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)

include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS ArraySortTest DESTINATION bin/tests/)

add_test(ArraySortTest ArraySortTest)
