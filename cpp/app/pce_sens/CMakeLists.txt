project (UQTk)

add_executable (pce_sens pce_sens.cpp)

target_link_libraries (pce_sens uqtkpce  )
target_link_libraries (pce_sens uqtkquad )
target_link_libraries (pce_sens uqtktools)
target_link_libraries (pce_sens uqtkarray)

target_link_libraries (pce_sens depdsfmt )
target_link_libraries (pce_sens depcvode )
target_link_libraries (pce_sens depnvec  )
target_link_libraries (pce_sens depslatec)
target_link_libraries (pce_sens deplapack)
target_link_libraries (pce_sens depblas  )
target_link_libraries (pce_sens depfigtree )
target_link_libraries (pce_sens depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_sens gfortran stdc++)
  else()
    target_link_libraries (pce_sens ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_sens ifcore)
  else()
    target_link_libraries (pce_sens ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_sens gfortran stdc++)
  else()
    target_link_libraries (pce_sens ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS pce_sens DESTINATION bin)

