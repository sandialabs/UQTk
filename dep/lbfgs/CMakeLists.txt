
enable_language(Fortran)

add_library(deplbfgs lbfgsDR.c lbfgs_routines.f)
include_directories("../../cpp/lib/include")

INSTALL(TARGETS deplbfgs DESTINATION lib)

SET(lbfgs_HEADERS
  lbfgs_routines.h
  )

INSTALL(FILES ${lbfgs_HEADERS} DESTINATION include/dep)
