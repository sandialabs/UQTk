
SET(kle_HEADERS
  kle.h
  )

add_library(uqtkkle kle.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )

include_directories (../../../dep/lapack)

# Install the library
INSTALL(TARGETS uqtkkle DESTINATION lib)

# Install the header files
INSTALL(FILES ${kle_HEADERS} DESTINATION include/uqtk)
