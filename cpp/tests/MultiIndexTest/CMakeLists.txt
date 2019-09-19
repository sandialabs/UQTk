project (UQTk)

add_executable (MultiIndexTest main.cpp)

target_link_libraries (MultiIndexTest uqtk)

target_link_libraries (MultiIndexTest depdsfmt   )
target_link_libraries (MultiIndexTest sundials_cvode)
target_link_libraries (MultiIndexTest sundials_nvecserial)
target_link_libraries (MultiIndexTest sundials_sunlinsoldense)
target_link_libraries (MultiIndexTest sundials_sunmatrixdense)
target_link_libraries (MultiIndexTest depslatec  )
target_link_libraries (MultiIndexTest m blas ${BLAS_LIBRARIES})
target_link_libraries (MultiIndexTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (MultiIndexTest deplbfgs   )
target_link_libraries (MultiIndexTest depfigtree )
target_link_libraries (MultiIndexTest depann     )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (MultiIndexTest gfortran stdc++)
  else()
    target_link_libraries (MultiIndexTest ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (MultiIndexTest ifcore)
  else()
    target_link_libraries (MultiIndexTest ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (MultiIndexTest gfortran stdc++)
  else()
    target_link_libraries (MultiIndexTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()


include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/mcmc  )


include_directories(../../../dep/dsfmt)
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/figtree)

INSTALL(TARGETS MultiIndexTest DESTINATION bin/tests/)

add_test(MultiIndexTest MultiIndexTest)
