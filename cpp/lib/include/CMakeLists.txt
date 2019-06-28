
SET(misc_HEADERS
  error_handlers.h
  ftndefs.h
  gen_defs.h
  uqtkconfig.h
  )

# Install the header files
INSTALL(FILES ${misc_HEADERS} DESTINATION include/uqtk)
