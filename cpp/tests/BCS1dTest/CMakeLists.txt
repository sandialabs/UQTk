project (UQTk)

add_executable (BCS1dTest main.cpp)

target_link_libraries (BCS1dTest uqtk )

target_link_libraries (BCS1dTest depdsfmt )
target_link_libraries (BCS1dTest sundials_cvode)
target_link_libraries (BCS1dTest sundials_nvecserial)
target_link_libraries (BCS1dTest sundials_sunlinsoldense)
target_link_libraries (BCS1dTest sundials_sunmatrixdense)
target_link_libraries (BCS1dTest depslatec)
target_link_libraries (BCS1dTest deplbfgs  )
target_link_libraries (BCS1dTest depfigtree  )
target_link_libraries (BCS1dTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (BCS1dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (BCS1dTest ifcore ifport)
  else()
    target_link_libraries (BCS1dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (BCS1dTest ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (BCS1dTest gfortran stdc++)
  else()
    target_link_libraries (BCS1dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries(BCS1dTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries(BCS1dTest m blas ${BLAS_LIBRARIES})

include_directories(../../lib/pce    )
include_directories(../../lib/bcs    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/mcmc  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()

INSTALL(TARGETS BCS1dTest DESTINATION bin/tests/)

add_test(BCS1dTest BCS1dTest)
