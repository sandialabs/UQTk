project (UQTk)

INSTALL(FILES ${copy_FILES} DESTINATION cpp/tests/MCMC2dTest/)

add_executable (MCMC2dTest main.cpp)

target_link_libraries (MCMC2dTest uqtk  )

target_link_libraries (MCMC2dTest depdsfmt )
target_link_libraries (MCMC2dTest depcvode )
target_link_libraries (MCMC2dTest depnvec  )
target_link_libraries (MCMC2dTest depslatec)
target_link_libraries (MCMC2dTest deplapack)
target_link_libraries (MCMC2dTest depblas  )
target_link_libraries (MCMC2dTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (MCMC2dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (MCMC2dTest ifcore ifport)
  else()
    target_link_libraries (MCMC2dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (MCMC2dTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (MCMC2dTest gfortran stdc++)
  else()
    target_link_libraries (MCMC2dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS MCMC2dTest DESTINATION bin/tests/)

add_test(MCMC2dTest MCMC2dTest)
