
add_executable (gen_mi gen_mi.cpp)

target_link_libraries (gen_mi uqtkpce  )
target_link_libraries (gen_mi uqtkarray)
target_link_libraries (gen_mi uqtktools) 

target_link_libraries (gen_mi depdsfmt ) 
target_link_libraries (gen_mi deplapack)
target_link_libraries (gen_mi depblas  )
target_link_libraries (gen_mi depfigtree )
target_link_libraries (gen_mi depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (gen_mi gfortran stdc++)
  else()
    target_link_libraries (gen_mi ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (gen_mi ifcore)
  else()
    target_link_libraries (gen_mi ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (gen_mi gfortran stdc++)
  else()
    target_link_libraries (gen_mi ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/include)
include_directories(../../lib/array  )
include_directories(../../lib/tools  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)


INSTALL(TARGETS gen_mi DESTINATION bin)

