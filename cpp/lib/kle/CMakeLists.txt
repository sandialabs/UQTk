SET(kle_HEADERS
  kle.h
  )

add_library(uqtkkle kle.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )

#find and include lapack
target_link_libraries(uqtkkle m lapack ${LAPACK_LIBRARIES})

# Install the library
INSTALL(TARGETS uqtkkle DESTINATION lib)

# Install the header files
INSTALL(FILES ${kle_HEADERS} DESTINATION include/uqtk)
