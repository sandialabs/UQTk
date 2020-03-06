project(UQTk)

SET(quad_HEADERS
  quad.h
  )

add_library(uqtkquad quad.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )  

# Install the library
INSTALL(TARGETS uqtkquad DESTINATION lib)

# Install the header files
INSTALL(FILES ${quad_HEADERS} DESTINATION include/uqtk)
