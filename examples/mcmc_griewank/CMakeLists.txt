project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  driver.cpp
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  cluster.sh
  SpectralClustering.py
  )

SET(copy_EXE_FILES
    mcmc_griewank.x
    calcLPGriewank.x
    griewankLL.x
   )

add_executable (griewankLL.x model.cpp)
add_executable (calcLPGriewank.x lprior.cpp)
target_link_libraries (calcLPGriewank.x uqtktmcmc depdsfmt)

# mcmc_griewank.x
add_executable (mcmc_griewank.x driver.cpp)
target_link_libraries (mcmc_griewank.x uqtktmcmc depdsfmt depslatec deplapack depblas)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (mcmc_griewank.x  gfortran stdc++)
  else()
    target_link_libraries (mcmc_griewank.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (mcmc_griewank.x    ifcore ifport)
  else()
    target_link_libraries (mcmc_griewank.x  ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (mcmc_griewank.x gfortran stdc++)
  else()
    target_link_libraries (mcmc_griewank.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/mcmc_griewank)
include_directories(../../cpp/lib/include )
include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_griewank
)
INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_griewank
)
