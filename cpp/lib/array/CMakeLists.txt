project(UQTk)

SET(array_HEADERS
  Array1D.h
  Array2D.h
  Array3D.h
  arrayio.h
  arraytools.h
  )

add_library(uqtkarray arrayio.cpp arraytools.cpp)

include_directories (../include)

#find and include blas and lapack
target_link_libraries(uqtkarray blas lapack)

# Install the library
INSTALL(TARGETS uqtkarray DESTINATION lib)

# Install the header files
INSTALL(FILES ${array_HEADERS} DESTINATION include/uqtk)
