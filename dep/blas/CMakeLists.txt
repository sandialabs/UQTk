
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(blas ${slsrc})

INSTALL(TARGETS blas DESTINATION lib)

# Install the header files
SET(blas_HEADERS
  depblas.h
  )

INSTALL(FILES ${blas_HEADERS} DESTINATION include/dep)
