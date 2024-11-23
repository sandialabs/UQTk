project (UQTk)

add_executable (BCS2dTest main.cpp)

# SET(copy_FILES
#   mindex_2d.dat
#   ck_2d.dat
#   )
# INSTALL(FILES ${copy_FILES} DESTINATION bin/tests/)

target_link_libraries (BCS2dTest uqtk )

target_link_libraries (BCS2dTest depdsfmt )
target_link_libraries (BCS2dTest sundials_core)
target_link_libraries (BCS2dTest sundials_cvode)
target_link_libraries (BCS2dTest sundials_nvecserial)
target_link_libraries (BCS2dTest sundials_sunlinsoldense)
target_link_libraries (BCS2dTest sundials_sunmatrixdense)
target_link_libraries (BCS2dTest depslatec)
target_link_libraries (BCS2dTest deplbfgs  )
target_link_libraries (BCS2dTest depfigtree  )
target_link_libraries (BCS2dTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (BCS2dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (BCS2dTest ifcore ifport)
  else()
    target_link_libraries (BCS2dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (BCS2dTest ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (BCS2dTest gfortran stdc++)
  else()
    target_link_libraries (BCS2dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (BCS2dTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (BCS2dTest m blas ${BLAS_LIBRARIES})

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

INSTALL(TARGETS BCS2dTest DESTINATION bin/tests/)

add_test(BCS2dTest BCS2dTest)
