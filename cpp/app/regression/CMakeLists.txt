project (UQTk)

add_executable (regression regression.cpp)

target_link_libraries (regression uqtklreg )
target_link_libraries (regression uqtkpce )
target_link_libraries (regression uqtkbcs  )
target_link_libraries (regression uqtkquad )
target_link_libraries (regression uqtktools)
target_link_libraries (regression uqtkarray)

target_link_libraries (regression depdsfmt )
target_link_libraries (regression deplbfgs )
target_link_libraries (regression sundials_cvode)
target_link_libraries (regression sundials_nvecserial)
target_link_libraries (regression sundials_sunlinsoldense)
target_link_libraries (regression sundials_sunmatrixdense)
target_link_libraries (regression depslatec)
target_link_libraries (regression depfigtree )
target_link_libraries (regression depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (regression gfortran stdc++)
  else()
    target_link_libraries (regression ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (regression ifcore)
    target_link_libraries (regression ifport)
  else()
    target_link_libraries (regression ${IntelLibPath}/libifcore.a)
    target_link_libraries (regression ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (regression gfortran stdc++)
  else()
    target_link_libraries (regression ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (regression m lapack ${LAPACK_LIBRARIES})
target_link_libraries (regression m blas ${BLAS_LIBRARIES})

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/bcs    )
include_directories(../../lib/lreg   )


include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/figtree)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories (../../../dep/sundials/include )
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
endif()

INSTALL(TARGETS regression DESTINATION bin)
