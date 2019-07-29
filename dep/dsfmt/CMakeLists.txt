
add_library(depdsfmt dSFMT.c dSFMT_h.c dsfmt_add.c)

INSTALL(TARGETS depdsfmt DESTINATION lib)

# Install the header files
SET(dsfmt_HEADERS
  dSFMT-params.h
  dSFMT-params11213.h
  dSFMT-params1279.h
  dSFMT-params132049.h
  dSFMT-params19937.h
  dSFMT-params216091.h
  dSFMT-params2203.h
  dSFMT-params4253.h
  dSFMT-params44497.h
  dSFMT-params521.h
  dSFMT-params86243.h
  dSFMT.h
  dsfmt_add.h
  )

INSTALL(FILES ${dsfmt_HEADERS} DESTINATION include/dep)

