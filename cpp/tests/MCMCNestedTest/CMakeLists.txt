project (UQTk)

INSTALL(FILES ${copy_FILES} DESTINATION cpp/tests/MCMCNestedTest/)

add_executable (MCMCNestedTest main.cpp)

target_link_libraries (MCMCNestedTest uqtk  )

target_link_libraries (MCMCNestedTest depdsfmt )
target_link_libraries (MCMCNestedTest depcvode )
target_link_libraries (MCMCNestedTest depnvec  )
target_link_libraries (MCMCNestedTest depslatec)
target_link_libraries (MCMCNestedTest deplapack)
target_link_libraries (MCMCNestedTest depblas  )
target_link_libraries (MCMCNestedTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (MCMCNestedTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (MCMCNestedTest ifcore ifport)
  else()
    target_link_libraries (MCMCNestedTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (MCMCNestedTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (MCMCNestedTest gfortran stdc++)
  else()
    target_link_libraries (MCMCNestedTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS MCMCNestedTest DESTINATION bin/tests/)

add_test(MCMCNestedTest MCMCNestedTest)
