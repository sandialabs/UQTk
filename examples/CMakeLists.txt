project (UQTk)

add_subdirectory (ops)
add_subdirectory (kle_ex1)
add_subdirectory (line_infer)
add_subdirectory (pce_bcs)
add_subdirectory (surf_rxn)
add_subdirectory (uqpc)
add_subdirectory (iuq)
add_subdirectory (sensMC)
add_subdirectory (sampling_qmc)
add_subdirectory (mcmc_rb)
add_subdirectory (mcmc_pol)
add_subdirectory (mcmc_griewank)
add_subdirectory (mcmc_manifold)
add_subdirectory (mcmc_hdf5)
add_subdirectory (num_integ)
add_subdirectory (fwd_prop)
add_subdirectory (window)
add_subdirectory (polynomial)
add_subdirectory (d_spring_series)


# Add muq example only with the gnu compilers
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND EXISTS ${PATH2MUQ})
   # message(STATUS "Adding MUQ<->UQTk example")
   #add_subdirectory (muq)
endif()

# Copy over README file too 
INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/)
