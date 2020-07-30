project(UQTk)

enable_language(Fortran)

SET(mmala_HEADERS
  mmala.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtkmmala ${slsrc})
 
include_directories (../../../dep/dsfmt )
include_directories (../tmcmc )
include_directories (../mcmc )
include_directories (../mala )
include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/figtree)
include_directories (../../../dep/lbfgs )
include_directories (../../../dep/slatec)


# Install the header files
INSTALL(FILES ${uqtkmmala_HEADERS} DESTINATION include/uqtk)
