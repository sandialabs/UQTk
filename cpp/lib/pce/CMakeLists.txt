project(UQTk)

SET(pce_HEADERS
  PCSet.h
  PCBasis.h
  )

add_library(uqtkpce PCBasis.cpp PCSet.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )

include_directories (../../../dep/slatec)
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/cvode-2.7.0/include)  
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

# Install the library
INSTALL(TARGETS uqtkpce DESTINATION lib)

# Install the header files
INSTALL(FILES ${pce_HEADERS} DESTINATION include/uqtk)
