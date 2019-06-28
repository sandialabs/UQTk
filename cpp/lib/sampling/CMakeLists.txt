project(UQTk)

SET(sampling_HEADERS
  sampling.hpp
   )
 
add_library(uqtksampling lhs.cpp qmc.cpp)

include_directories (../include)
include_directories (../array)
include_directories (../tools)

include_directories (../../../dep/dsfmt)
include_directories (../../../dep/figtree)

# Install the library
INSTALL(TARGETS uqtksampling DESTINATION lib)

# Install the header files
INSTALL(FILES ${sampling_HEADERS} DESTINATION include/uqtk)
