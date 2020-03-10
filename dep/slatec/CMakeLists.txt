
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(depslatec ${slsrc})

INSTALL(TARGETS depslatec DESTINATION lib)

SET(slatec_HEADERS
  depslatec.h
  )

INSTALL(FILES ${slatec_HEADERS} DESTINATION include/dep)
