project (UQTk)

add_executable (pce_quad pce_quad.cpp)

target_link_libraries (pce_quad uqtkpce  )
target_link_libraries (pce_quad uqtkquad )
target_link_libraries (pce_quad uqtktools)
target_link_libraries (pce_quad uqtkarray)

target_link_libraries (pce_quad depdsfmt )
target_link_libraries (pce_quad depcvode )
target_link_libraries (pce_quad depnvec  )
target_link_libraries (pce_quad depslatec)
target_link_libraries (pce_quad deplapack)
target_link_libraries (pce_quad depblas  )
target_link_libraries (pce_quad depfigtree )
target_link_libraries (pce_quad depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_quad gfortran stdc++)
  else()
    target_link_libraries (pce_quad ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_quad ifcore)
  else()
    target_link_libraries (pce_quad ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_quad gfortran stdc++)
  else()
    target_link_libraries (pce_quad ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS pce_quad DESTINATION bin)

