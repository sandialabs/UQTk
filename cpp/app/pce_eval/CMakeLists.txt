project (UQTk)

add_executable (pce_eval pce_eval.cpp)

target_link_libraries (pce_eval uqtkpce  )
target_link_libraries (pce_eval uqtkquad )
target_link_libraries (pce_eval uqtktools)
target_link_libraries (pce_eval uqtkarray)

target_link_libraries (pce_eval depdsfmt )
target_link_libraries (pce_eval sundials_cvode)
target_link_libraries (pce_eval sundials_nvecserial)
target_link_libraries (pce_eval sundials_sunlinsoldense)
target_link_libraries (pce_eval sundials_sunmatrixdense)
target_link_libraries (pce_eval depslatec)
target_link_libraries (pce_eval depfigtree )
target_link_libraries (pce_eval depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_eval gfortran stdc++)
  else()
    target_link_libraries (pce_eval ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_eval ifcore)
  else()
    target_link_libraries (pce_eval ${IntelLibPath}/libifcore.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_eval gfortran stdc++)
  else()
    target_link_libraries (pce_eval ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (pce_eval m lapack ${LAPACK_LIBRARIES})
target_link_libraries (pce_eval m blas ${BLAS_LIBRARIES})

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )


include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories (../../../dep/sundials/include )
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
endif()
INSTALL(TARGETS pce_eval DESTINATION bin)
