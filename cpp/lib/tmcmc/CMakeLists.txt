project(UQTk)

enable_language(Fortran)

SET(tmcmc_HEADERS
  tmcmc.h
  )

FILE(GLOB slsrc "*.cpp")

add_library(uqtktmcmc ${slsrc})
 
include_directories (../../../dep/dsfmt )

# Install the header files
INSTALL(FILES ${uqtktmcmc_HEADERS} DESTINATION include/uqtk)

