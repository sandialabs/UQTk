project(UQTk)

enable_language(Fortran)

SET(mala_HEADERS
  mala.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtkmala ${slsrc})
 
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
INSTALL(FILES ${uqtkmala_HEADERS} DESTINATION include/uqtk)
