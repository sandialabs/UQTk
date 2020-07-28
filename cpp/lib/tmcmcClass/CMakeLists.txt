project(UQTk)

enable_language(Fortran)

SET(tmcmcClass_HEADERS
  tmcmcClass.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtktmcmcClass ${slsrc})
 
include_directories (../../../dep/dsfmt )
include_directories (../tmcmc )
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
INSTALL(FILES ${uqtktmcmcClass_HEADERS} DESTINATION include/uqtk)

