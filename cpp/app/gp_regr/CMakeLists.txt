project (UQTk)

add_executable (gp_regr gp_regr.cpp)

target_link_libraries (gp_regr uqtkgproc )
target_link_libraries (gp_regr uqtkpce  )
target_link_libraries (gp_regr uqtkbcs  )
target_link_libraries (gp_regr uqtkquad )
target_link_libraries (gp_regr uqtktools)
target_link_libraries (gp_regr uqtkarray)

target_link_libraries (gp_regr depdsfmt )
target_link_libraries (gp_regr deplbfgs )
target_link_libraries (gp_regr depcvode )
target_link_libraries (gp_regr depnvec  )
target_link_libraries (gp_regr depslatec)
target_link_libraries (gp_regr deplapack)
target_link_libraries (gp_regr depblas  )
target_link_libraries (gp_regr depfigtree )
target_link_libraries (gp_regr depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (gp_regr gfortran stdc++)
  else()
    target_link_libraries (gp_regr ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (gp_regr ifcore)
    target_link_libraries (gp_regr ifport)
  else()
    target_link_libraries (gp_regr ${IntelLibPath}/libifcore.a)
    target_link_libraries (gp_regr ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (gp_regr gfortran stdc++)
  else()
    target_link_libraries (gp_regr ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/bcs    )
include_directories(../../lib/gproc  )


include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS gp_regr DESTINATION bin)

