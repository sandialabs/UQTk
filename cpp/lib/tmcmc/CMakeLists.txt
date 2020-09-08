project(UQTk)

enable_language(Fortran)

SET(tmcmc_HEADERS
  tmcmc.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtktmcmc ${slsrc})
 
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

# Install the header files
INSTALL(FILES ${uqtktmcmc_HEADERS} DESTINATION include/uqtk)

