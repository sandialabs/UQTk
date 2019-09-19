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


target_link_libraries(uqtklowrank m blas ${BLAS_LIBRARIES})
target_link_libraries(uqtklowrank m lapack ${LAPACK_LIBRARIES})
include_directories (../../../dep/lbfgs)
include_directories (../../../dep/dsfmt)
include_directories (../../../dep/figtree)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()

# Install the library
INSTALL(TARGETS uqtklowrank DESTINATION lib)

# Install the header files
INSTALL(FILES ${lowrank_HEADERS} DESTINATION include/uqtk)
