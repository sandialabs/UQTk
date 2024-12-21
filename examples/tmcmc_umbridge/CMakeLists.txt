project (UQTk)

SET(copy_FILES
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  tmcmc_moveIntermediateFiles.sh
  tmcmc_prior_samples.dat
  tmcmc_umbridge.py
  minimal-umbridge-model.py
  README
  )

# link_directories(${CMAKE_SUNDIALS_DIR}/lib)

add_executable (umbridge_model.x umbridge_model.cpp)
target_link_libraries (umbridge_model.x pthread)

add_executable (tmcmc_umbridge.x tmcmc_umbridge.cpp)
target_link_libraries (tmcmc_umbridge.x pthread)
target_link_libraries (tmcmc_umbridge.x uqtk  )
target_link_libraries (tmcmc_umbridge.x depdsfmt )
target_link_libraries (tmcmc_umbridge.x depslatec)
target_link_libraries (tmcmc_umbridge.x deplbfgs  )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (tmcmc_umbridge.x gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (tmcmc_umbridge.x ifcore ifport)
  else()
    target_link_libraries (tmcmc_umbridge.x ${IntelLibPath}/libifcore.a)
    target_link_libraries (tmcmc_umbridge.x ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (tmcmc_umbridge.x gfortran stdc++)
  else()
    target_link_libraries (tmcmc_umbridge.x ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (tmcmc_umbridge.x m lapack)
target_link_libraries (tmcmc_umbridge.x m blas)

include_directories(../../examples/tmcmc_umbridge)

include_directories(../../cpp/lib/pce    )
include_directories(../../cpp/lib/array  )
include_directories(../../cpp/lib/include)
include_directories(../../cpp/lib/quad   )
include_directories(../../cpp/lib/tools  )
include_directories(../../cpp/lib/mcmc  )
include_directories(../../cpp/lib/tmcmc  )

include_directories(../../dep/dsfmt)
include_directories(../../dep/slatec)
include_directories(../../dep/lbfgs)


INSTALL(TARGETS tmcmc_umbridge.x DESTINATION examples/tmcmc_umbridge)
INSTALL(TARGETS umbridge_model.x DESTINATION examples/tmcmc_umbridge)

INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/tmcmc_umbridge)