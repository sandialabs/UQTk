SET(tmcmc_HEADERS
  tmcmc.h
  )

add_library(uqtktmcmc tmcmc.cpp)
 
include_directories (../../../dep/dsfmt )
include_directories (../mcmc )
include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../tmcmc  )
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/figtree)
include_directories (../../../dep/lbfgs )
include_directories (../../../dep/slatec)

include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()

# Install the library
INSTALL(TARGETS uqtktmcmc DESTINATION lib)

# Install the header files
INSTALL(FILES ${tmcmc_HEADERS} DESTINATION include/uqtk)

