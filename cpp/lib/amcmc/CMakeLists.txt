project(UQTk)

enable_language(Fortran)

SET(amcmc_HEADERS
  amcmc.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtkamcmc ${slsrc})
 
include_directories (../../../dep/dsfmt )
include_directories (../tmcmc )
include_directories (../mcmc )
include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/figtree)
include_directories (../../../dep/lbfgs )
include_directories (../../../dep/slatec)


# Install the header files
INSTALL(FILES ${uqtkamcmc_HEADERS} DESTINATION include/uqtk)
