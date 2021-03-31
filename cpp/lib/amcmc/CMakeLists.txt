SET(amcmc_HEADERS
  amcmc.h
  )

add_library(uqtkamcmc amcmc.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../mcmc  )

include_directories (../../../dep/slatec)
#target_link_libraries(uqtkmcmc m lapack ${LAPACK_LIBRARIES})
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/figtree)
include_directories (../../../dep/lbfgs )

include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
  include_directories (../../../dep/sundials/include )
endif()



# Install the library
INSTALL(TARGETS uqtkamcmc DESTINATION lib)

# Install the header files
INSTALL(FILES ${amcmc_HEADERS} DESTINATION include/uqtk)
