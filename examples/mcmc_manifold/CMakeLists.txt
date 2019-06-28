project (UQTk)

SET(copy_SCRIPT_FILES
  PyKDE.py
  run.sh
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  move.sh
  setup.dat
  rngs.dat
  )

SET(copy_EXE_FILES
    mcmc_manifold.x
    calcLPmanifold.x
    manifoldLL.x
   )

add_executable (manifoldLL.x model.cpp)
add_executable (calcLPmanifold.x lprior.cpp)
target_link_libraries (calcLPmanifold.x uqtktmcmc depdsfmt)

# mcmc_manifold.x
add_executable (mcmc_manifold.x driver.cpp)
target_link_libraries (mcmc_manifold.x uqtktmcmc depdsfmt depslatec deplapack depblas)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (mcmc_manifold.x  gfortran stdc++)
  else()
    target_link_libraries (mcmc_manifold.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (mcmc_manifold.x    ifcore ifport)
  else()
    target_link_libraries (mcmc_manifold.x  ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (mcmc_manifold.x gfortran stdc++)
  else()
    target_link_libraries (mcmc_manifold.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/mcmc_manifold)
include_directories(../../cpp/lib/include )
include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_manifold
)
INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_manifold
)
