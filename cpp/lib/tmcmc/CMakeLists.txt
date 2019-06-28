
SET(tmcmc_HEADERS
  tmcmc.h
  )

add_library(uqtktmcmc tmcmc.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )

include_directories (../../../dep/lapack)
include_directories (../../../dep/dsfmt )
#if (NOT "${HDF5LibPath}" STREQUAL "")
  include_directories(${HDF5LibPath}/include)
#endif()


include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

# Install the library
INSTALL(TARGETS uqtktmcmc DESTINATION lib)

# Install the header files
INSTALL(FILES ${tmcmc_HEADERS} DESTINATION include/uqtk)
