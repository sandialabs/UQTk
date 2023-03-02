project (UQTk)

add_executable (PCELogTest main.cpp)

target_link_libraries (PCELogTest uqtk )

target_link_libraries (PCELogTest depdsfmt )
target_link_libraries (PCELogTest sundials_cvode)
target_link_libraries (PCELogTest sundials_nvecserial)
target_link_libraries (PCELogTest sundials_sunlinsoldense)
target_link_libraries (PCELogTest sundials_sunmatrixdense)
target_link_libraries (PCELogTest depslatec)
target_link_libraries (PCELogTest deplbfgs  )
target_link_libraries (PCELogTest depfigtree  )
target_link_libraries (PCELogTest depann  )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (PCELogTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (PCELogTest ifcore ifport)
  else()
    target_link_libraries (PCELogTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (PCELogTest ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (PCELogTest gfortran stdc++)
  else()
    target_link_libraries (PCELogTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (PCELogTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (PCELogTest m blas ${BLAS_LIBRARIES})

include_directories(../../lib/pce    )
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

INSTALL(TARGETS PCELogTest DESTINATION bin/tests/)

add_test(PCELogTest PCELogTest)
