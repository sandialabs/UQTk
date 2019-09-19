project (UQTk)

SET(copy_SCRIPT_FILES
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  tmcmc_moveIntermediateFiles.sh
  )

SET(copy_DATA_FILES
  tmcmc_prior_samples.dat
  )

SET(copy_EXE_FILES
    model.x
   )

file(COPY ${copy_SCRIPT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${copy_DATA_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable (model.x model.cpp)
add_executable (TMCMC2dTest main.cpp)

target_link_libraries (TMCMC2dTest uqtk  )

target_link_libraries (TMCMC2dTest depdsfmt )
target_link_libraries (TMCMC2dTest sundials_cvode)
target_link_libraries (TMCMC2dTest sundials_nvecserial)
target_link_libraries (TMCMC2dTest sundials_sunlinsoldense)
target_link_libraries (TMCMC2dTest sundials_sunmatrixdense)
target_link_libraries (TMCMC2dTest depslatec)
target_link_libraries (TMCMC2dTest m lapack blas)
target_link_libraries (TMCMC2dTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (TMCMC2dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (TMCMC2dTest ifcore ifport)
  else()
    target_link_libraries (TMCMC2dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (TMCMC2dTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (TMCMC2dTest gfortran stdc++)
  else()
    target_link_libraries (TMCMC2dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS TMCMC2dTest DESTINATION bin/tests/)

add_test(TMCMC2dTest TMCMC2dTest)
