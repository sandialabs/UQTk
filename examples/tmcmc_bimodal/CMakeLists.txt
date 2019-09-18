project (UQTk)

SET(copy_FILES
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  tmcmc_moveIntermediateFiles.sh
  tmcmc_prior_samples.dat
  tmcmc_bimodal.py
  )

add_executable (bimodal.x bimodal.cpp)
add_executable (tmcmc_bimodal.x tmcmc_bimodal.cpp)

target_link_libraries (tmcmc_bimodal.x uqtk  )
target_link_libraries (tmcmc_bimodal.x depdsfmt )
target_link_libraries (tmcmc_bimodal.x depcvode )
target_link_libraries (tmcmc_bimodal.x depnvec  )
target_link_libraries (tmcmc_bimodal.x depslatec)
target_link_libraries (tmcmc_bimodal.x m lapack blas)
target_link_libraries (tmcmc_bimodal.x deplbfgs  )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (tmcmc_bimodal.x gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (tmcmc_bimodal.x ifcore ifport)
  else()
    target_link_libraries (tmcmc_bimodal.x ${IntelLibPath}/libifcore.a)
    target_link_libraries (tmcmc_bimodal.x ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (tmcmc_bimodal.x gfortran stdc++)
  else()
    target_link_libraries (tmcmc_bimodal.x ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/tmcmc_bimodal)

include_directories(../../cpp/lib/pce    )
include_directories(../../cpp/lib/array  )
include_directories(../../cpp/lib/include)
include_directories(../../cpp/lib/quad   )
include_directories(../../cpp/lib/tools  )
include_directories(../../cpp/lib/mcmc  )

include_directories(../../dep/dsfmt)
include_directories(../../dep/slatec)
include_directories(../../dep/lbfgs)
include_directories(../../dep/cvode-2.7.0/include)

INSTALL(TARGETS tmcmc_bimodal.x DESTINATION examples/tmcmc_bimodal)
INSTALL(TARGETS bimodal.x DESTINATION examples/tmcmc_bimodal)

INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/tmcmc_bimodal)
