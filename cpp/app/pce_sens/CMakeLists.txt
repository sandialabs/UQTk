project (UQTk)

add_executable (pce_sens pce_sens.cpp)

target_link_libraries (pce_sens uqtkpce  )
target_link_libraries (pce_sens uqtkquad )
target_link_libraries (pce_sens uqtktools)
target_link_libraries (pce_sens uqtkarray)

target_link_libraries (pce_sens depdsfmt )
target_link_libraries (pce_sens sundials_cvode)
target_link_libraries (pce_sens sundials_nvecserial)
target_link_libraries (pce_sens sundials_sunlinsoldense)
target_link_libraries (pce_sens sundials_sunmatrixdense)
target_link_libraries (pce_sens depslatec)
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
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_sens gfortran stdc++)
  else()
    target_link_libraries (pce_sens ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (pce_sens m lapack ${LAPACK_LIBRARIES})
target_link_libraries (pce_sens m blas ${BLAS_LIBRARIES})

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

INSTALL(TARGETS pce_sens DESTINATION bin)
