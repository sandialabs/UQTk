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
find_package(BLAS)
find_package(LAPACK)
target_link_libraries(uqtkarray blas lapack)
#target_link_libraries(uqtkarray m blas ${BLAS_LIBRARIES})
#target_link_libraries(uqtkarray m lapack ${LAPACK_LIBRARIES})

# Install the library
INSTALL(TARGETS uqtkarray DESTINATION lib)

# Install the header files
INSTALL(FILES ${array_HEADERS} DESTINATION include/uqtk)
