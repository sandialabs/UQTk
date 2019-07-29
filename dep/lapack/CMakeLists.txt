
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(deplapack ${slsrc})

INSTALL(TARGETS deplapack DESTINATION lib)

SET(lapack_HEADERS
  deplapack.h
  )

INSTALL(FILES ${lapack_HEADERS} DESTINATION include/dep)
