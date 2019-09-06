project(UQTk)

SET(lowrank_HEADERS
  lowrank.h
  )

add_library(uqtklowrank lowrank.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../pce    )
include_directories (../bcs    )


include_directories (../../../dep/lapack)
include_directories (../../../dep/blas)
include_directories (../../../dep/lbfgs)
include_directories (../../../dep/dsfmt)
include_directories (../../../dep/figtree)
include_directories (../../../dep/sundials/include )
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (${CMAKE_SUNDIALS_DIR}/include)

# Install the library
INSTALL(TARGETS uqtklowrank DESTINATION lib)

# Install the header files
INSTALL(FILES ${lowrank_HEADERS} DESTINATION include/uqtk)
