project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  move.sh
  SpectralClustering.py
  hdf5clean.sh
  )

SET(copy_EXE_FILES
    modelLL.x
    modelLP.x
    mcmc_hdf5.x
   )

add_executable (modelLL.x model.cpp)
add_executable (modelLP.x lprior.cpp)
target_link_libraries (modelLP.x uqtktmcmc depdsfmt)


# mcmc_hdf5.x
add_executable (mcmc_hdf5.x driver.cpp)
target_link_libraries (mcmc_hdf5.x uqtktmcmc depdsfmt depslatec deplapack depblas)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (mcmc_hdf5.x  gfortran stdc++)
  else()
    target_link_libraries (mcmc_hdf5.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (mcmc_hdf5.x    ifcore ifport)
  else()
    target_link_libraries (mcmc_hdf5.x  ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (mcmc_hdf5.x gfortran stdc++)
  else()
    target_link_libraries (mcmc_hdf5.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

if (NOT "${HDF5LibPath}" STREQUAL "")
  target_link_libraries (mcmc_hdf5.x  ${HDF5LibPath}/lib/libhdf5.a ${HDF5LibPath}/lib/libhdf5_cpp.a ${HDF5LibPath}/lib/libhdf5_hl.a ${HDF5LibPath}/lib/libhdf5_hl_cpp.a)
endif()

include_directories(../../examples/mcmc_hdf5)
include_directories(../../cpp/lib/include )
include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )

if (NOT "${HDF5LibPath}" STREQUAL "")
  include_directories(${HDF5LibPath}/include)
endif()

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_hdf5
)
INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/mcmc_hdf5
)
