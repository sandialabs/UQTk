project (UQTk)

add_executable (ArrayDelColumn main.cpp)

target_link_libraries (ArrayDelColumn uqtk )

target_link_libraries (ArrayDelColumn depdsfmt )
target_link_libraries (ArrayDelColumn depcvode )
target_link_libraries (ArrayDelColumn depnvec  )
target_link_libraries (ArrayDelColumn depslatec)
target_link_libraries (ArrayDelColumn deplapack)
target_link_libraries (ArrayDelColumn depblas  )
target_link_libraries (ArrayDelColumn deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (ArrayDelColumn gfortran expat stdc++) 
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (ArrayDelColumn ifcore)
  else()
    target_link_libraries (ArrayDelColumn ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (ArrayDelColumn gfortran stdc++)
  else()
    target_link_libraries (ArrayDelColumn ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS ArrayDelColumn DESTINATION bin/tests/)

add_test(ArrayDelColumn ArrayDelColumn)
