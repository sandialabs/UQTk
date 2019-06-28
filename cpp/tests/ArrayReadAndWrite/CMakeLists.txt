project (UQTk)

add_executable (ArrayReadAndWrite main.cpp)

target_link_libraries (ArrayReadAndWrite uqtk  )

target_link_libraries (ArrayReadAndWrite depdsfmt )
target_link_libraries (ArrayReadAndWrite depcvode )
target_link_libraries (ArrayReadAndWrite depnvec  )
target_link_libraries (ArrayReadAndWrite depslatec)
target_link_libraries (ArrayReadAndWrite deplapack)
target_link_libraries (ArrayReadAndWrite depblas  )
target_link_libraries (ArrayReadAndWrite deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (ArrayReadAndWrite gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (ArrayReadAndWrite ifcore)
  else()
    target_link_libraries (ArrayReadAndWrite ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (ArrayReadAndWrite gfortran expat stdc++)
  else()
    target_link_libraries (ArrayReadAndWrite ${ClangLibPath}/libgfortran.dylib expat ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS ArrayReadAndWrite DESTINATION bin/tests/)

add_test(ArrayReadAndWrite ArrayReadAndWrite)
