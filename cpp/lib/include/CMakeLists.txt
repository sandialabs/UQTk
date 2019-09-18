
SET(misc_HEADERS
  error_handlers.h
  ftndefs.h
  gen_defs.h
  uqtkconfig.h
  depblas.h
  deplapack.h
  )

# Install the header files
INSTALL(FILES ${misc_HEADERS} DESTINATION include/uqtk)
