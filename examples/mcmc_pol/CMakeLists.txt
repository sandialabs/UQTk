project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  getEv.sh
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  model.py
  yd_n11_s0.01.dat
  )

SET(copy_EXE_FILES
    mcmc_pol.x
    polLP.x
  )

# mcmc_pol.x
add_executable (mcmc_pol.x driver.cpp)
add_executable (polLP.x lprior.cpp)
target_link_libraries (mcmc_pol.x uqtktmcmc depdsfmt depslatec deplapack depblas)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (mcmc_pol.x  gfortran stdc++)
  else()
    target_link_libraries (mcmc_pol.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (mcmc_pol.x    ifcore ifport)
  else()
    target_link_libraries (mcmc_pol.x  ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (mcmc_pol.x gfortran stdc++)
  else()
    target_link_libraries (mcmc_pol.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/mcmc_pol)
include_directories(../../cpp/lib/include )
include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_pol
)
INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_pol
)

