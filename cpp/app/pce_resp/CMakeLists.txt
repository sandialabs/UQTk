project (UQTk)

add_executable (pce_resp pce_resp.cpp)

target_link_libraries (pce_resp uqtkpce  )
target_link_libraries (pce_resp uqtkquad )
target_link_libraries (pce_resp uqtktools)
target_link_libraries (pce_resp uqtkarray)

target_link_libraries (pce_resp depdsfmt )
target_link_libraries (pce_resp depcvode )
target_link_libraries (pce_resp depnvec  )
target_link_libraries (pce_resp depslatec)
target_link_libraries (pce_resp deplapack)
target_link_libraries (pce_resp depblas  )
target_link_libraries (pce_resp depfigtree )
target_link_libraries (pce_resp depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_resp gfortran stdc++)
  else()
    target_link_libraries (pce_resp ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_resp ifcore)
  else()
    target_link_libraries (pce_resp ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_resp gfortran stdc++)
  else()
    target_link_libraries (pce_resp ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS pce_resp DESTINATION bin)

