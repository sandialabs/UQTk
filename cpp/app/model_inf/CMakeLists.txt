project (UQTk)

add_executable (model_inf model_inf.cpp)

target_link_libraries (model_inf uqtkinfer)
target_link_libraries (model_inf uqtkmcmc )
target_link_libraries (model_inf uqtkpce  )
target_link_libraries (model_inf uqtkquad )
target_link_libraries (model_inf uqtktools)
target_link_libraries (model_inf uqtkarray)

target_link_libraries (model_inf depdsfmt )
target_link_libraries (model_inf deplbfgs )
target_link_libraries (model_inf depcvode )
target_link_libraries (model_inf depnvec  )
target_link_libraries (model_inf depslatec)
target_link_libraries (model_inf deplapack)
target_link_libraries (model_inf depblas  )
target_link_libraries (model_inf depfigtree )
target_link_libraries (model_inf depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (model_inf gfortran stdc++)
  else()
    target_link_libraries (model_inf ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (model_inf ifcore ifport)
  else()
    target_link_libraries (model_inf ${IntelLibPath}/libifcore.a)
    target_link_libraries (model_inf ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (model_inf gfortran stdc++)
  else()
    target_link_libraries (model_inf ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce     )
include_directories(../../lib/array   )
include_directories(../../lib/include )
include_directories(../../lib/quad    )
include_directories(../../lib/tools   )
include_directories(../../lib/mcmc    )
include_directories(../../lib/infer   )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")


INSTALL(TARGETS model_inf DESTINATION bin)

