project (UQTk)

INSTALL(FILES ${copy_FILES} DESTINATION cpp/tests/MCMCRandomTest/)

add_executable (MCMCRandomTest main.cpp)

target_link_libraries (MCMCRandomTest uqtk  )

target_link_libraries (MCMCRandomTest depdsfmt )
target_link_libraries (MCMCRandomTest sundials_cvode)
target_link_libraries (MCMCRandomTest sundials_nvecserial)
target_link_libraries (MCMCRandomTest sundials_sunlinsoldense)
target_link_libraries (MCMCRandomTest sundials_sunmatrixdense)
target_link_libraries (MCMCRandomTest depslatec)
target_link_libraries (MCMCRandomTest m blas ${BLAS_LIBRARIES})
target_link_libraries (MCMCRandomTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (MCMCRandomTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (MCMCRandomTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (MCMCRandomTest ifcore ifport)
  else()
    target_link_libraries (MCMCRandomTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (MCMCRandomTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (MCMCRandomTest gfortran stdc++)
  else()
    target_link_libraries (MCMCRandomTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS MCMCRandomTest DESTINATION bin/tests/)

add_test(MCMCRandomTest MCMCRandomTest)
