
add_executable (generate_quad generate_quad.cpp)

target_link_libraries (generate_quad uqtkquad )
target_link_libraries (generate_quad uqtktools)
target_link_libraries (generate_quad uqtkarray)

target_link_libraries (generate_quad depdsfmt )
target_link_libraries (generate_quad deplapack)
target_link_libraries (generate_quad depblas  )
target_link_libraries (generate_quad depfigtree )
target_link_libraries (generate_quad depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (generate_quad gfortran stdc++)
  else()
    target_link_libraries (generate_quad ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (generate_quad ifcore)
  else()
    target_link_libraries (generate_quad ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (generate_quad gfortran stdc++)
  else()
    target_link_libraries (generate_quad ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/include)
include_directories(../../lib/tools  )
include_directories(../../lib/quad   )
include_directories(../../lib/array  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)


INSTALL(TARGETS generate_quad DESTINATION bin)

