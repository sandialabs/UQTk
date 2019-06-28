project(UQTk)

SET(bcs_HEADERS
  bcs.h
   )
 
add_library(uqtkbcs bcs.cpp)

include_directories (../include)
include_directories (../array)
include_directories (../tools)
include_directories (../quad)

include_directories (../../../dep/slatec)
include_directories (../../../dep/lapack) 
include_directories (../../../dep/dsfmt)
include_directories (../../../dep/figtree)
include_directories (../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

# Install the library
INSTALL(TARGETS uqtkbcs DESTINATION lib)

# Install the header files
INSTALL(FILES ${bcs_HEADERS} DESTINATION include/uqtk)
