project(UQTk)

enable_language(Fortran)

SET(dfi_HEADERS
  dfi.h
  )

add_library(uqtkdfi dfi.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../mcmc  )
include_directories (../quad  )
include_directories (../sampling  )
include_directories (../tools  )
include_directories (../../../dep/figtree  )
include_directories (../pce  )

include_directories (../../../dep/dsfmt )
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()


# Install the library
INSTALL(TARGETS uqtkdfi DESTINATION lib)

# Install the header files
INSTALL(FILES ${dfi_HEADERS} DESTINATION include/uqtk)
