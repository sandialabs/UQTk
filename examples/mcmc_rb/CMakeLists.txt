project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  driver.cpp
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  rngs.dat
  )

SET(copy_EXE_FILES
    mcmc_rb.x
    calcLP.x
   )


add_executable (calcLP.x lprior.cpp)
target_link_libraries (calcLP.x uqtktmcmc depdsfmt)

# mcmc_rb.x
add_executable (mcmc_rb.x driver.cpp)
target_link_libraries (mcmc_rb.x uqtktmcmc depdsfmt depslatec deplapack depblas)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (mcmc_rb.x  gfortran stdc++)
  else()
    target_link_libraries (mcmc_rb.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (mcmc_rb.x    ifcore ifport)
  else()
    target_link_libraries (mcmc_rb.x  ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (mcmc_rb.x gfortran stdc++)
  else()
    target_link_libraries (mcmc_rb.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/mcmc_rb)
include_directories(../../cpp/lib/include )
include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_rb
)
INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_rb
)

