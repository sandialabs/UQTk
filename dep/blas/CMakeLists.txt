
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(depblas ${slsrc})

INSTALL(TARGETS depblas DESTINATION lib)

# Install the header files
SET(blas_HEADERS
  depblas.h
  )

INSTALL(FILES ${blas_HEADERS} DESTINATION include/dep)

