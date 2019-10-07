project (UQTk)

add_executable (Array1DMiscTest main.cpp)

target_link_libraries (Array1DMiscTest uqtk )

target_link_libraries (Array1DMiscTest depdsfmt )
target_link_libraries (Array1DMiscTest sundials_cvode)
target_link_libraries (Array1DMiscTest sundials_nvecserial)
target_link_libraries (Array1DMiscTest sundials_sunlinsoldense)
target_link_libraries (Array1DMiscTest sundials_sunmatrixdense)
target_link_libraries (Array1DMiscTest depslatec)
target_link_libraries (Array1DMiscTest m blas ${BLAS_LIBRARIES})
target_link_libraries (Array1DMiscTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (Array1DMiscTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (Array1DMiscTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (Array1DMiscTest ifcore ifport)
  else()
    target_link_libraries (Array1DMiscTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (Array1DMiscTest ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (Array1DMiscTest gfortran stdc++)
  else()
    target_link_libraries (Array1DMiscTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS Array1DMiscTest DESTINATION bin/tests/)

add_test(Array1DMiscTest Array1DMiscTest)
