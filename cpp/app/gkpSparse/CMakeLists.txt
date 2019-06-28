enable_language(Fortran)
enable_language(CXX)

add_executable (gkpSparse gkpSparse.cpp gkpclib.cpp gkpflib.f)

target_link_libraries (gkpSparse uqtktools)
target_link_libraries (gkpSparse uqtkarray)

target_link_libraries (gkpSparse depdsfmt )
target_link_libraries (gkpSparse deplapack)
target_link_libraries (gkpSparse depblas  )
target_link_libraries (gkpSparse depfigtree )
target_link_libraries (gkpSparse depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (gkpSparse gfortran stdc++)
  else()
    target_link_libraries (gkpSparse ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (gkpSparse ifcore ifport)
  else()
    target_link_libraries (gkpSparse ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (gkpSparse gfortran stdc++)
  else()
    target_link_libraries (gkpSparse ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/include)
include_directories(../../lib/tools  )
include_directories(../../lib/array  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)

INSTALL(TARGETS gkpSparse DESTINATION bin)

