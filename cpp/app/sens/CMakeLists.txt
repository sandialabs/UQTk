
add_executable (sens sens.cpp)
add_executable (trdSpls trdSpls.cpp)

target_link_libraries (trdSpls uqtkmcmc )
target_link_libraries (trdSpls uqtkpce  )
target_link_libraries (trdSpls uqtkquad )
target_link_libraries (trdSpls uqtktools)
target_link_libraries (trdSpls uqtkarray)

target_link_libraries (trdSpls depdsfmt )
target_link_libraries (trdSpls deplbfgs )
target_link_libraries (trdSpls depcvode )
target_link_libraries (trdSpls depnvec  )
target_link_libraries (trdSpls depslatec)
target_link_libraries (trdSpls deplapack)
target_link_libraries (trdSpls depblas  )
target_link_libraries (trdSpls depfigtree )
target_link_libraries (trdSpls depann   )


target_link_libraries (sens uqtkarray)

target_link_libraries (sens depdsfmt )
target_link_libraries (sens deplapack)
target_link_libraries (sens depblas  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (sens gfortran stdc++)
    target_link_libraries (trdSpls gfortran stdc++)
  else()
    target_link_libraries (sens ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
    target_link_libraries (trdSpls ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (sens ifcore)
    target_link_libraries (trdSpls ifcore)
  else()
    target_link_libraries (sens ${IntelLibPath}/libifcore.a)
    target_link_libraries (trdSpls ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (sens gfortran stdc++)
    target_link_libraries (trdSpls gfortran stdc++)
  else()
    target_link_libraries (sens ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
    target_link_libraries (trdSpls ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/include)
include_directories(../../lib/tools  )
include_directories(../../lib/quad   )
include_directories(../../lib/array  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lapack)
include_directories(../../../dep/blas)
include_directories(../../../dep/figtree)

INSTALL(TARGETS sens DESTINATION bin)
INSTALL(TARGETS trdSpls DESTINATION bin)

