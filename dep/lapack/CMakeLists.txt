
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(lapack ${slsrc})

INSTALL(TARGETS lapack DESTINATION lib)

SET(lapack_HEADERS
  deplapack.h
  )

INSTALL(FILES ${lapack_HEADERS} DESTINATION include/dep)
